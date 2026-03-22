from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import networkx as nx
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

from graynet.writers import (
    RAW_EDGE_FILE_NAME,
    ROI_STATS_FILE_NAME,
    SUMMARY_EDGE_FILE_NAME,
    read_run_metadata,
)


def _add_node_positions(graph: nx.Graph, metadata: dict) -> None:
    centroids = metadata.get("centroids", {})
    for node, coords in centroids.items():
        graph.add_node(node)
        graph.nodes[node]["x"] = float(coords[0])
        graph.nodes[node]["y"] = float(coords[1])
        graph.nodes[node]["z"] = float(coords[2])


def _group_graph_rows(rows, group_keys):
    grouped = defaultdict(list)
    for row in rows:
        grouped[tuple(row[key] for key in group_keys)].append(row)
    return grouped


def export_graphml(run_dir, out_dir=None) -> Path:
    run_dir = Path(run_dir).resolve()
    export_dir = Path(out_dir).resolve() if out_dir is not None else run_dir / "exports" / "graphml"
    export_dir.mkdir(exist_ok=True, parents=True)
    metadata = read_run_metadata(run_dir)

    raw_path = run_dir / RAW_EDGE_FILE_NAME
    if raw_path.exists():
        raw_rows = pq.read_table(raw_path).to_pylist()
        for (subject_id, base_feature, weight_method), rows in _group_graph_rows(
            raw_rows, ["subject_id", "base_feature", "weight_method"]
        ).items():
            graph = nx.Graph()
            _add_node_positions(graph, metadata)
            for row in rows:
                graph.add_edge(row["u"], row["v"], weight=float(row["weight"]))
            subject_dir = export_dir / subject_id
            subject_dir.mkdir(exist_ok=True, parents=True)
            out_path = subject_dir / f"{base_feature}__{weight_method}.graphml"
            nx.write_graphml(graph, out_path, encoding="utf-8")

    summary_path = run_dir / SUMMARY_EDGE_FILE_NAME
    if summary_path.exists():
        summary_rows = pq.read_table(summary_path).to_pylist()
        for (subject_id, weight_method, summary_stat), rows in _group_graph_rows(
            summary_rows, ["subject_id", "weight_method", "summary_stat"]
        ).items():
            graph = nx.Graph()
            _add_node_positions(graph, metadata)
            for row in rows:
                graph.add_edge(row["u"], row["v"], weight=float(row["weight"]))
            subject_dir = export_dir / subject_id
            subject_dir.mkdir(exist_ok=True, parents=True)
            out_path = subject_dir / f"{weight_method}__{summary_stat}__summary.graphml"
            nx.write_graphml(graph, out_path, encoding="utf-8")

    return export_dir


def export_csv(run_dir, out_dir=None) -> Path:
    run_dir = Path(run_dir).resolve()
    export_dir = Path(out_dir).resolve() if out_dir is not None else run_dir / "exports" / "csv"
    export_dir.mkdir(exist_ok=True, parents=True)

    for file_name in (RAW_EDGE_FILE_NAME, SUMMARY_EDGE_FILE_NAME, ROI_STATS_FILE_NAME):
        source = run_dir / file_name
        if source.exists():
            table = pq.read_table(source)
            pacsv.write_csv(table, export_dir / file_name.replace(".parquet", ".csv"))

    return export_dir
