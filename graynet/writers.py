from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from graynet.domain import RunConfig

RAW_EDGE_FILE_NAME = "edges_raw.parquet"
SUMMARY_EDGE_FILE_NAME = "edges_summary.parquet"
ROI_STATS_FILE_NAME = "roi_stats.parquet"
METADATA_FILE_NAME = "run_metadata.json"

RAW_EDGE_SCHEMA = pa.schema([
    ("subject_id", pa.string()),
    ("base_feature", pa.string()),
    ("weight_method", pa.string()),
    ("u", pa.string()),
    ("v", pa.string()),
    ("weight", pa.float64()),
])

SUMMARY_EDGE_SCHEMA = pa.schema([
    ("subject_id", pa.string()),
    ("weight_method", pa.string()),
    ("summary_stat", pa.string()),
    ("u", pa.string()),
    ("v", pa.string()),
    ("weight", pa.float64()),
])

ROI_STATS_SCHEMA = pa.schema([
    ("subject_id", pa.string()),
    ("base_feature", pa.string()),
    ("stat_name", pa.string()),
    ("roi", pa.string()),
    ("value", pa.float64()),
])


def _slugify(value: Any) -> str:
    text = str(value).strip().lower()
    keep = []
    for char in text:
        if char.isalnum():
            keep.append(char)
        elif char in {"-", "_", "."}:
            keep.append(char)
        else:
            keep.append("-")
    collapsed = "".join(keep).strip("-")
    while "--" in collapsed:
        collapsed = collapsed.replace("--", "-")
    return collapsed or "none"


def _joined(values: tuple[str, ...] | list[str]) -> str:
    return "+".join(_slugify(value) for value in values) if values else "none"


def build_run_dir(config: RunConfig) -> Path | None:
    if config.out_dir is None:
        return None

    parts = [
        _slugify(config.mode),
        f"atlas-{_slugify(config.atlas_name)}",
        f"features-{_joined(list(config.base_features))}",
    ]
    if config.weight_methods:
        parts.append(f"weights-{_joined(list(config.weight_methods))}")
    if config.summary_stat_names:
        parts.append(f"summaries-{_joined(list(config.summary_stat_names))}")
    if config.roi_stat_names:
        parts.append(f"stats-{_joined(list(config.roi_stat_names))}")
    parts.append(f"smth-{_slugify(config.smoothing_param)}")
    parts.append(f"node-{_slugify(config.node_size)}")

    run_dir = config.out_dir / "__".join(parts)
    run_dir.mkdir(exist_ok=True, parents=True)
    return run_dir


class ParquetBatchWriter:
    def __init__(self, path: Path, schema: pa.Schema):
        self.path = Path(path)
        self.schema = schema
        self._writer: pq.ParquetWriter | None = None

    def write_rows(self, rows: list[dict[str, Any]]) -> None:
        if not rows:
            return

        table = pa.Table.from_pylist(rows, schema=self.schema)
        if self._writer is None:
            self._writer = pq.ParquetWriter(str(self.path), self.schema)
        self._writer.write_table(table)

    def close(self) -> None:
        if self._writer is not None:
            self._writer.close()
            self._writer = None

    def __enter__(self) -> "ParquetBatchWriter":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()


def write_run_metadata(run_dir: Path, metadata: dict[str, Any]) -> Path:
    metadata_path = run_dir / METADATA_FILE_NAME
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    return metadata_path


def read_run_metadata(run_dir: Path) -> dict[str, Any]:
    metadata_path = Path(run_dir) / METADATA_FILE_NAME
    return json.loads(metadata_path.read_text())
