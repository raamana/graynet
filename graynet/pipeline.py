from __future__ import annotations

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path
from typing import Any

import hiwenet
import networkx as nx
import numpy as np

from graynet import config_graynet as cfg
from graynet.atlas import atlas_identifier, resolve_atlas
from graynet.domain import AtlasInfo, RunConfig, SubjectBatchResult
from graynet.utils import (
    calc_roi_statistics,
    check_edge_range_dict,
    check_features,
    check_num_bins,
    check_num_procs,
    check_params_multiedge,
    check_params_single_edge,
    check_stat_methods,
    check_subjects,
    check_weight_params,
    check_weights,
    import_features,
    mask_background_roi,
    warn_nan,
)
from graynet.writers import (
    METADATA_FILE_NAME,
    RAW_EDGE_FILE_NAME,
    RAW_EDGE_SCHEMA,
    ROI_STATS_FILE_NAME,
    ROI_STATS_SCHEMA,
    SUMMARY_EDGE_FILE_NAME,
    SUMMARY_EDGE_SCHEMA,
    ParquetBatchWriter,
    build_run_dir,
    write_run_metadata,
)


def _node_label(node) -> str:
    if isinstance(node, str):
        return node

    try:
        value = float(node)
    except (TypeError, ValueError):
        return str(node)

    if value.is_integer():
        return str(int(value))

    return str(value)


def _vector_from_rows(rows: list[dict[str, Any]], node_order: tuple[str, ...]) -> np.ndarray:
    node_index = {node: idx for idx, node in enumerate(node_order)}
    sorted_rows = sorted(rows, key=lambda row: (node_index[row["u"]], node_index[row["v"]]))
    return np.array([row["weight"] for row in sorted_rows], dtype=float)


def _rows_from_graph(
    graph: nx.Graph,
    subject_id: str,
    base_feature: str,
    weight_method: str,
    node_order: tuple[str, ...],
) -> list[dict[str, Any]]:
    node_index = {node: idx for idx, node in enumerate(node_order)}
    rows = []
    for u, v, attrs in graph.edges(data=True):
        u_label = _node_label(u)
        v_label = _node_label(v)
        if node_index[u_label] > node_index[v_label]:
            u_label, v_label = v_label, u_label
        rows.append(
            {
                "subject_id": subject_id,
                "base_feature": base_feature,
                "weight_method": weight_method,
                "u": u_label,
                "v": v_label,
                "weight": float(attrs["weight"]),
            }
        )
    return rows


def _load_masked_feature(
    input_dir: Path,
    subject_id: str,
    base_feature: str,
    smoothing_param,
    atlas_spec,
    roi_labels,
    ignore_label,
):
    features = import_features(
        input_dir,
        [subject_id],
        base_feature,
        fwhm=smoothing_param,
        atlas=atlas_spec,
    )
    return mask_background_roi(features[subject_id], roi_labels, ignore_label)


def _run_single_subject_edges(
    subject_id: str,
    input_dir: Path,
    base_feature: str,
    smoothing_param,
    atlas_spec,
    atlas_info: AtlasInfo,
    weight_methods: tuple[str, ...],
    num_bins: int,
    edge_range,
    return_results: bool,
) -> SubjectBatchResult:
    raw_rows: list[dict[str, Any]] = []
    payload: dict[Any, Any] = {}
    data, rois = _load_masked_feature(
        input_dir,
        subject_id,
        base_feature,
        smoothing_param,
        atlas_spec,
        atlas_info.roi_labels,
        atlas_info.ignore_label,
    )

    for weight_method in weight_methods:
        graph = hiwenet.extract(
            data,
            rois,
            weight_method=weight_method,
            num_bins=num_bins,
            edge_range=edge_range,
            return_networkx_graph=True,
        )
        rows = _rows_from_graph(graph, subject_id, base_feature, weight_method, atlas_info.node_labels)
        vector = _vector_from_rows(rows, atlas_info.node_labels)
        warn_nan(vector)
        raw_rows.extend(rows)
        if return_results:
            payload[(weight_method, subject_id)] = vector

    return SubjectBatchResult(raw_rows, [], [], payload)


def _run_single_subject_roi_stats(
    subject_id: str,
    input_dir: Path,
    base_feature: str,
    smoothing_param,
    atlas_spec,
    atlas_info: AtlasInfo,
    stat_funcs: tuple[Any, ...],
    stat_names: tuple[str, ...],
    return_results: bool,
) -> SubjectBatchResult:
    data, rois = _load_masked_feature(
        input_dir,
        subject_id,
        base_feature,
        smoothing_param,
        atlas_spec,
        atlas_info.roi_labels,
        atlas_info.ignore_label,
    )

    roi_stat_rows: list[dict[str, Any]] = []
    payload: dict[Any, Any] = {}
    unique_rois = atlas_info.roi_values
    roi_lookup = np.asarray(unique_rois)

    for stat_func, stat_name in zip(stat_funcs, stat_names):
        roi_values = calc_roi_statistics(data, rois, roi_lookup, stat_func)
        for roi, value in zip(atlas_info.node_labels, roi_values):
            roi_stat_rows.append(
                {
                    "subject_id": subject_id,
                    "base_feature": base_feature,
                    "stat_name": stat_name,
                    "roi": roi,
                    "value": float(value),
                }
            )
        if return_results:
            key = subject_id if len(stat_names) == 1 else (stat_name, subject_id)
            payload[key] = np.asarray(roi_values, dtype=float)

    return SubjectBatchResult([], [], roi_stat_rows, payload)


def _run_single_subject_multiedge(
    subject_id: str,
    input_dir: Path,
    base_features: tuple[str, ...],
    smoothing_param,
    atlas_spec,
    atlas_info: AtlasInfo,
    weight_methods: tuple[str, ...],
    num_bins: int,
    edge_range_dict: dict[str, tuple[float, float] | None],
    summary_stats: tuple[Any, ...],
    summary_stat_names: tuple[str, ...],
    return_results: bool,
) -> SubjectBatchResult:
    feature_cache = {}
    for base_feature in base_features:
        feature_cache[base_feature] = _load_masked_feature(
            input_dir,
            subject_id,
            base_feature,
            smoothing_param,
            atlas_spec,
            atlas_info.roi_labels,
            atlas_info.ignore_label,
        )

    raw_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    payload: dict[Any, Any] = {}

    for weight_method in weight_methods:
        rows_for_weight: list[dict[str, Any]] = []
        for base_feature in base_features:
            data, rois = feature_cache[base_feature]
            graph = hiwenet.extract(
                data,
                rois,
                weight_method=weight_method,
                num_bins=num_bins,
                edge_range=edge_range_dict[base_feature],
                return_networkx_graph=True,
            )
            rows = _rows_from_graph(
                graph,
                subject_id,
                base_feature,
                weight_method,
                atlas_info.node_labels,
            )
            vector = _vector_from_rows(rows, atlas_info.node_labels)
            warn_nan(vector)
            rows_for_weight.extend(rows)
            raw_rows.extend(rows)
            if return_results:
                payload[(weight_method, base_feature, subject_id)] = vector

        grouped_weights: dict[tuple[str, str], list[float]] = defaultdict(list)
        for row in rows_for_weight:
            grouped_weights[(row["u"], row["v"])].append(row["weight"])

        for stat_func, stat_name in zip(summary_stats, summary_stat_names):
            for (u, v), weights in grouped_weights.items():
                summary_rows.append(
                    {
                        "subject_id": subject_id,
                        "weight_method": weight_method,
                        "summary_stat": stat_name,
                        "u": u,
                        "v": v,
                        "weight": float(stat_func(np.asarray(weights, dtype=float))),
                    }
                )

    return SubjectBatchResult(raw_rows, summary_rows, [], payload)


def _collect_results(
    worker,
    subject_ids: tuple[str, ...],
    num_procs: int,
) -> list[SubjectBatchResult]:
    if num_procs <= 1 or len(subject_ids) <= 1:
        return [worker(subject_id) for subject_id in subject_ids]

    with ProcessPoolExecutor(max_workers=num_procs) as pool:
        return list(pool.map(worker, subject_ids))


def _single_edge_worker(
    subject_id: str,
    input_dir: Path,
    base_feature: str,
    smoothing_param,
    atlas_spec,
    atlas_info: AtlasInfo,
    weight_methods: tuple[str, ...],
    num_bins: int,
    edge_range,
    return_results: bool,
) -> SubjectBatchResult:
    return _run_single_subject_edges(
        subject_id,
        input_dir,
        base_feature,
        smoothing_param,
        atlas_spec,
        atlas_info,
        weight_methods,
        num_bins,
        edge_range,
        return_results,
    )


def _roi_stats_worker(
    subject_id: str,
    input_dir: Path,
    base_feature: str,
    smoothing_param,
    atlas_spec,
    atlas_info: AtlasInfo,
    stat_funcs: tuple[Any, ...],
    stat_names: tuple[str, ...],
    return_results: bool,
) -> SubjectBatchResult:
    return _run_single_subject_roi_stats(
        subject_id,
        input_dir,
        base_feature,
        smoothing_param,
        atlas_spec,
        atlas_info,
        stat_funcs,
        stat_names,
        return_results,
    )


def _multiedge_worker(
    subject_id: str,
    input_dir: Path,
    base_features: tuple[str, ...],
    smoothing_param,
    atlas_spec,
    atlas_info: AtlasInfo,
    weight_methods: tuple[str, ...],
    num_bins: int,
    edge_range_dict: dict[str, tuple[float, float] | None],
    summary_stats: tuple[Any, ...],
    summary_stat_names: tuple[str, ...],
    return_results: bool,
) -> SubjectBatchResult:
    return _run_single_subject_multiedge(
        subject_id,
        input_dir,
        base_features,
        smoothing_param,
        atlas_spec,
        atlas_info,
        weight_methods,
        num_bins,
        edge_range_dict,
        summary_stats,
        summary_stat_names,
        return_results,
    )


def _write_metadata(run_dir: Path, config: RunConfig, atlas_info: AtlasInfo) -> None:
    metadata = {
        "mode": config.mode,
        "atlas_name": atlas_info.atlas_name,
        "atlas_spec": atlas_identifier(config.atlas_spec, atlas_info.atlas_name),
        "base_features": list(config.base_features),
        "weight_methods": list(config.weight_methods),
        "summary_stat_names": list(config.summary_stat_names),
        "roi_stat_names": list(config.roi_stat_names),
        "subject_ids": list(config.subject_ids),
        "node_labels": list(atlas_info.node_labels),
        "centroids": {node: list(coords) for node, coords in atlas_info.centroids.items()},
        "smoothing_param": config.smoothing_param,
        "node_size": config.node_size,
        "num_bins": config.num_bins,
        "edge_range": list(config.edge_range) if config.edge_range is not None else None,
        "edge_range_dict": config.edge_range_dict,
        "files": {
            "metadata": METADATA_FILE_NAME,
            "raw_edges": RAW_EDGE_FILE_NAME,
            "summary_edges": SUMMARY_EDGE_FILE_NAME if config.summary_stat_names else None,
            "roi_stats": ROI_STATS_FILE_NAME if config.roi_stat_names else None,
        },
    }
    write_run_metadata(run_dir, metadata)


def _run_with_writers(
    config: RunConfig,
    atlas_info: AtlasInfo,
    worker,
) -> tuple[dict[Any, Any], Path | None]:
    run_dir = build_run_dir(config)
    if run_dir is not None:
        _write_metadata(run_dir, config, atlas_info)

    raw_writer = (
        ParquetBatchWriter(run_dir / RAW_EDGE_FILE_NAME, RAW_EDGE_SCHEMA) if run_dir else None
    )
    summary_writer = (
        ParquetBatchWriter(run_dir / SUMMARY_EDGE_FILE_NAME, SUMMARY_EDGE_SCHEMA)
        if run_dir and config.summary_stat_names
        else None
    )
    roi_writer = (
        ParquetBatchWriter(run_dir / ROI_STATS_FILE_NAME, ROI_STATS_SCHEMA)
        if run_dir and config.roi_stat_names
        else None
    )

    payload: dict[Any, Any] = {}
    writers = [writer for writer in (raw_writer, summary_writer, roi_writer) if writer]
    try:
        for batch in _collect_results(worker, config.subject_ids, config.num_procs):
            if raw_writer:
                raw_writer.write_rows(batch.raw_edge_rows)
            if summary_writer:
                summary_writer.write_rows(batch.summary_edge_rows)
            if roi_writer:
                roi_writer.write_rows(batch.roi_stat_rows)
            payload.update(batch.return_payload)
    finally:
        for writer in writers:
            writer.close()

    return payload, run_dir


def run_edges(
    subject_id_list,
    input_dir,
    base_feature=cfg.default_feature_single_edge,
    weight_method_list=cfg.default_weight_method,
    num_bins=cfg.default_num_bins,
    edge_range=cfg.default_edge_range,
    atlas=cfg.default_atlas,
    smoothing_param=cfg.default_smoothing_param,
    node_size=cfg.default_node_size,
    out_dir=None,
    return_results=False,
    num_procs=cfg.default_num_procs,
):
    feature_name = check_features(base_feature)[0]
    check_params_single_edge(feature_name, input_dir, atlas, smoothing_param, node_size, out_dir, return_results)
    atlas_info = resolve_atlas(feature_name, atlas, node_size)
    subject_ids, _, _, _ = check_subjects(subject_id_list)
    weight_methods, _, _, _ = check_weights(weight_method_list)
    num_bins, edge_range = check_weight_params(num_bins, edge_range)
    num_procs = check_num_procs(num_procs)

    config = RunConfig(
        mode="edges",
        input_dir=Path(input_dir),
        out_dir=Path(out_dir) if out_dir is not None else None,
        subject_ids=tuple(str(subject) for subject in subject_ids),
        base_features=(feature_name,),
        weight_methods=tuple(weight_methods),
        atlas_spec=atlas_info.atlas_spec,
        atlas_name=atlas_info.atlas_name,
        smoothing_param=smoothing_param,
        node_size=node_size,
        num_bins=num_bins,
        edge_range=edge_range,
        return_results=return_results,
        num_procs=num_procs,
    )

    worker = partial(
        _single_edge_worker,
        input_dir=config.input_dir,
        base_feature=feature_name,
        smoothing_param=smoothing_param,
        atlas_spec=atlas_info.atlas_spec,
        atlas_info=atlas_info,
        weight_methods=config.weight_methods,
        num_bins=num_bins,
        edge_range=edge_range,
        return_results=return_results,
    )
    payload, run_dir = _run_with_writers(config, atlas_info, worker)
    return payload if return_results else run_dir


def run_roi_stats(
    subject_id_list,
    input_dir,
    base_feature=cfg.default_feature_single_edge,
    chosen_roi_stats=cfg.default_roi_statistic,
    atlas=cfg.default_atlas,
    smoothing_param=cfg.default_smoothing_param,
    node_size=cfg.default_node_size,
    out_dir=None,
    return_results=False,
    num_procs=cfg.default_num_procs,
):
    feature_name = check_features(base_feature)[0]
    check_params_single_edge(feature_name, input_dir, atlas, smoothing_param, node_size, out_dir, return_results)
    stat_funcs, stat_names, _, _, _ = check_stat_methods(chosen_roi_stats)
    atlas_info = resolve_atlas(feature_name, atlas, node_size)
    subject_ids, _, _, _ = check_subjects(subject_id_list)
    num_procs = check_num_procs(num_procs)

    config = RunConfig(
        mode="roi-stats",
        input_dir=Path(input_dir),
        out_dir=Path(out_dir) if out_dir is not None else None,
        subject_ids=tuple(str(subject) for subject in subject_ids),
        base_features=(feature_name,),
        weight_methods=(),
        atlas_spec=atlas_info.atlas_spec,
        atlas_name=atlas_info.atlas_name,
        smoothing_param=smoothing_param,
        node_size=node_size,
        num_bins=None,
        roi_stats=tuple(stat_funcs),
        roi_stat_names=tuple(stat_names),
        return_results=return_results,
        num_procs=num_procs,
    )

    worker = partial(
        _roi_stats_worker,
        input_dir=config.input_dir,
        base_feature=feature_name,
        smoothing_param=smoothing_param,
        atlas_spec=atlas_info.atlas_spec,
        atlas_info=atlas_info,
        stat_funcs=config.roi_stats,
        stat_names=config.roi_stat_names,
        return_results=return_results,
    )
    payload, run_dir = _run_with_writers(config, atlas_info, worker)
    return payload if return_results else run_dir


def run_multiedge(
    subject_id_list,
    input_dir,
    base_feature_list=cfg.default_features_multi_edge,
    weight_method_list=cfg.default_weight_method,
    summary_stats=cfg.multi_edge_summary_func_default,
    num_bins=cfg.default_num_bins,
    edge_range_dict=cfg.edge_range_predefined,
    atlas=cfg.default_atlas,
    smoothing_param=cfg.default_smoothing_param,
    node_size=cfg.default_node_size,
    out_dir=None,
    return_results=False,
    overwrite_results=False,
    num_procs=cfg.default_num_procs,
):
    for feature in base_feature_list:
        if feature in cfg.features_volumetric:
            raise NotImplementedError(
                "Multi-edge networks are not yet supported for volumetric features."
            )

    check_params_multiedge(base_feature_list, input_dir, atlas, smoothing_param, node_size, out_dir, return_results)
    features = tuple(check_features(base_feature_list))
    subject_ids, _, _, _ = check_subjects(subject_id_list)
    weight_methods, _, _, _ = check_weights(weight_method_list)
    stat_funcs, stat_names, _, _, _ = check_stat_methods(summary_stats)
    num_bins = check_num_bins(num_bins)
    given_ranges = dict(edge_range_dict) if edge_range_dict is not None else None
    edge_range_dict = check_edge_range_dict(given_ranges, features)
    num_procs = check_num_procs(num_procs)
    atlas_info = resolve_atlas(features[0], atlas, node_size)

    config = RunConfig(
        mode="multiedge",
        input_dir=Path(input_dir),
        out_dir=Path(out_dir) if out_dir is not None else None,
        subject_ids=tuple(str(subject) for subject in subject_ids),
        base_features=features,
        weight_methods=tuple(weight_methods),
        atlas_spec=atlas_info.atlas_spec,
        atlas_name=atlas_info.atlas_name,
        smoothing_param=smoothing_param,
        node_size=node_size,
        num_bins=num_bins,
        edge_range_dict=edge_range_dict,
        summary_stats=tuple(stat_funcs),
        summary_stat_names=tuple(stat_names),
        return_results=return_results,
        num_procs=num_procs,
    )

    worker = partial(
        _multiedge_worker,
        input_dir=config.input_dir,
        base_features=config.base_features,
        smoothing_param=smoothing_param,
        atlas_spec=atlas_info.atlas_spec,
        atlas_info=atlas_info,
        weight_methods=config.weight_methods,
        num_bins=num_bins,
        edge_range_dict=edge_range_dict,
        summary_stats=config.summary_stats,
        summary_stat_names=config.summary_stat_names,
        return_results=return_results,
    )
    payload, run_dir = _run_with_writers(config, atlas_info, worker)
    return payload if return_results else run_dir
