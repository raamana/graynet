from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator

import networkx as nx
import numpy as np
import pyarrow as pyarrow
import pyarrow.compute as pyarrow_compute
import pyarrow.parquet as pyarrow_parquet

from graynet.writers import (
    RAW_EDGE_FILE_NAME,
    ROI_STATS_FILE_NAME,
    SUMMARY_EDGE_FILE_NAME,
    read_run_metadata,
)


def _filter_table(table: pyarrow.Table, **filters) -> pyarrow.Table:
    mask = None
    for key, value in filters.items():
        if value is None:
            continue
        if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
            this_mask = pyarrow_compute.is_in(table[key], pyarrow.array(list(value)))
        else:
            this_mask = pyarrow_compute.equal(table[key], value)
        mask = this_mask if mask is None else pyarrow_compute.and_(mask, this_mask)

    return table if mask is None else table.filter(mask)


def _add_centroids(graph: nx.Graph, metadata: dict[str, Any]) -> None:
    for node, coords in metadata.get("centroids", {}).items():
        graph.add_node(node)
        graph.nodes[node]["x"] = float(coords[0])
        graph.nodes[node]["y"] = float(coords[1])
        graph.nodes[node]["z"] = float(coords[2])


@dataclass(frozen=True)
class EdgeData:
    table: pyarrow.Table
    metadata: dict[str, Any]

    def filter(
        self,
        subject_id: str | Sequence[str] | None = None,
        base_feature: str | None = None,
        weight_method: str | None = None,
        summary_stat: str | None = None,
    ) -> "EdgeData":
        return EdgeData(
            table=_filter_table(
                self.table,
                subject_id=subject_id,
                base_feature=base_feature,
                weight_method=weight_method,
                summary_stat=summary_stat,
            ),
            metadata=self.metadata,
        )

    def to_pandas(self):
        return self.table.to_pandas()

    def to_rows(self) -> list[dict[str, Any]]:
        return self.table.to_pylist()

    def stable_subject_ids(self) -> list[str]:
        """Subject IDs present in this table, in a defined (stabilized) order.

        Order is **not** arbitrary Parquet row order: IDs from run metadata
        ``subject_ids`` come first when they appear in the table (preserving
        that list's order), then any IDs found only in the table, in first-seen
        column order.

        .. note::

            This ordering is intentional for iteration and analysis, but it is
            **not** yet guaranteed as a stable public contract across graynet
            releases or if on-disk layout / metadata conventions change.
        """
        configured_order = list(self.metadata.get("subject_ids", []))
        present = set(self.table.column("subject_id").to_pylist())
        ordered = [subject_id for subject_id in configured_order if subject_id in present]

        configured_set = set(configured_order)
        extras = [
            subject_id
            for subject_id in dict.fromkeys(self.table.column("subject_id").to_pylist())
            if subject_id not in configured_set
        ]

        return ordered + extras

    def iter_subjects(self) -> Iterator[tuple[str, "EdgeData"]]:
        """Yield ``(subject_id, edge_rows_for_subject)`` using :meth:`stable_subject_ids` order."""
        for subject_id in self.stable_subject_ids():
            yield subject_id, self.filter(subject_id=subject_id)

    def to_ndarray(
        self,
        subject_ids: Sequence[str],
        base_feature: str | None = None,
        weight_method: str | None = None,
        summary_stat: str | None = None,
        fill_value: float = np.nan,
    ) -> np.ndarray:
        if isinstance(subject_ids, (str, bytes)):
            raise TypeError("subject_ids must be a sequence of subject identifiers, not a string")

        filtered = self.filter(
            base_feature=base_feature,
            weight_method=weight_method,
            summary_stat=summary_stat,
        )
        rows = filtered.to_rows()
        if not rows:
            return np.empty((len(subject_ids), 0), dtype=float)

        node_order = list(self.metadata.get("node_labels", []))
        node_index = {node: idx for idx, node in enumerate(node_order)}

        edge_keys = []
        seen_edges = set()
        for row in rows:
            edge_key = (row["u"], row["v"])
            if edge_key not in seen_edges:
                seen_edges.add(edge_key)
                edge_keys.append(edge_key)

        edge_keys.sort(key=lambda edge: (node_index.get(edge[0], float("inf")),
                                         node_index.get(edge[1], float("inf"))))
        edge_index = {edge: idx for idx, edge in enumerate(edge_keys)}
        subject_index = {subject_id: idx for idx, subject_id in enumerate(subject_ids)}

        matrix = np.full((len(subject_ids), len(edge_keys)), fill_value, dtype=float)
        for row in rows:
            sid = row["subject_id"]
            if sid not in subject_index:
                continue
            matrix[subject_index[sid], edge_index[(row["u"], row["v"])]] = float(row["weight"])

        return matrix


@dataclass(frozen=True)
class RoiStatsData:
    table: pyarrow.Table
    metadata: dict[str, Any]

    def filter(
        self,
        subject_id: str | None = None,
        base_feature: str | None = None,
        stat_name: str | None = None,
    ) -> "RoiStatsData":
        return RoiStatsData(
            table=_filter_table(
                self.table,
                subject_id=subject_id,
                base_feature=base_feature,
                stat_name=stat_name,
            ),
            metadata=self.metadata,
        )

    def to_pandas(self):
        return self.table.to_pandas()


@dataclass(frozen=True)
class RunData:
    run_dir: Path
    metadata: dict[str, Any]
    raw_edges: EdgeData | None = None
    summary_edges: EdgeData | None = None
    roi_stats: RoiStatsData | None = None

    def __iter__(self):
        yield self.raw_edges
        yield self.metadata


def load_run(run_dir) -> RunData:
    run_dir = Path(run_dir).resolve()
    metadata = read_run_metadata(run_dir)

    raw_edges = None
    raw_path = run_dir / RAW_EDGE_FILE_NAME
    if raw_path.exists():
        raw_edges = EdgeData(pyarrow_parquet.read_table(raw_path), metadata)

    summary_edges = None
    summary_path = run_dir / SUMMARY_EDGE_FILE_NAME
    if summary_path.exists():
        summary_edges = EdgeData(pyarrow_parquet.read_table(summary_path), metadata)

    roi_stats = None
    roi_path = run_dir / ROI_STATS_FILE_NAME
    if roi_path.exists():
        roi_stats = RoiStatsData(pyarrow_parquet.read_table(roi_path), metadata)

    return RunData(
        run_dir=run_dir,
        metadata=metadata,
        raw_edges=raw_edges,
        summary_edges=summary_edges,
        roi_stats=roi_stats,
    )


def get_edge_values(
    edge_data: EdgeData,
    subject_id: str | Sequence[str] | None = None,
    base_feature: str | None = None,
    weight_method: str | None = None,
    summary_stat: str | None = None,
) -> EdgeData:
    return edge_data.filter(
        subject_id=subject_id,
        base_feature=base_feature,
        weight_method=weight_method,
        summary_stat=summary_stat,
    )


def export_to_nx(edge_data: EdgeData, metadata: dict[str, Any] | None = None) -> nx.Graph:
    metadata = metadata or edge_data.metadata
    graph = nx.Graph()
    _add_centroids(graph, metadata)
    for row in edge_data.to_rows():
        graph.add_edge(row["u"], row["v"], weight=float(row["weight"]))
    return graph
