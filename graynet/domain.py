from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class AtlasInfo:
    atlas_spec: Any
    atlas_name: str
    roi_values: tuple[Any, ...]
    node_labels: tuple[str, ...]
    centroids: dict[str, tuple[float, float, float]]
    roi_labels: Any
    ignore_label: Any


@dataclass(frozen=True)
class RunConfig:
    mode: str
    input_dir: str | Path
    out_dir: str | Path | None
    subject_ids: Any
    atlas: Any
    smoothing_param: int | float | None
    node_size: int | None
    base_features: tuple[str, ...] = ()
    weight_methods: tuple[Any, ...] = ()
    atlas_name: str | None = None
    num_bins: int | None = None
    edge_range: tuple[float, float] | None = None
    edge_range_dict: dict[str, tuple[float, float] | None] | None = None
    summary_stats: tuple[Any, ...] = ()
    summary_stat_names: tuple[str, ...] = ()
    roi_stats: tuple[Any, ...] = ()
    roi_stat_names: tuple[str, ...] = ()
    return_results: bool = False
    num_procs: int = 1
    overwrite_results: bool = False


@dataclass(frozen=True)
class SubjectJob:
    subject_id: str
    config: RunConfig
    atlas_info: AtlasInfo


@dataclass
class SubjectResult:
    raw_edge_rows: list[dict[str, Any]]
    summary_edge_rows: list[dict[str, Any]]
    roi_stat_rows: list[dict[str, Any]]
    return_payload: dict[Any, Any]
