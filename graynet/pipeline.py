from __future__ import annotations

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path
from typing import Any

import hiwenet
import networkx as nx
import numpy as np

from graynet import config_graynet as cfg
from graynet.atlas import atlas_identifier, format_node_label, resolve_atlas
from graynet.domain import AtlasInfo, RunConfig, SubjectJob, SubjectResult
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


__all__ = ["run", "run_edges", "run_multiedge", "run_roi_stats"]


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
    rows: list[dict[str, Any]] = []
    for u, v, attrs in graph.edges(data=True):
        u_label = format_node_label(u)
        v_label = format_node_label(v)
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


def _load_feature(job: SubjectJob, base_feature: str | None = None):
    feature_name = base_feature or job.config.base_features[0]
    features = import_features(
        job.config.input_dir,
        [job.subject_id],
        feature_name,
        fwhm=job.config.smoothing_param,
        atlas=job.atlas_info.atlas_spec,
    )
    return mask_background_roi(
        features[job.subject_id],
        job.atlas_info.roi_labels,
        job.atlas_info.ignore_label,
    )


def _process_edges(job: SubjectJob) -> SubjectResult:
    config = job.config
    base_feature = config.base_features[0]
    data, rois = _load_feature(job)

    raw_rows: list[dict[str, Any]] = []
    payload: dict[Any, Any] = {}
    for weight_method in config.weight_methods:
        graph = hiwenet.extract(
            data,
            rois,
            weight_method=weight_method,
            num_bins=config.num_bins,
            edge_range=config.edge_range,
            return_networkx_graph=True,
        )
        rows = _rows_from_graph(
            graph,
            job.subject_id,
            base_feature,
            weight_method,
            job.atlas_info.node_labels,
        )
        vector = _vector_from_rows(rows, job.atlas_info.node_labels)
        warn_nan(vector)
        raw_rows.extend(rows)
        if config.return_results:
            payload[(weight_method, job.subject_id)] = vector

    return SubjectResult(raw_rows, [], [], payload)


def _process_roi_stats(job: SubjectJob) -> SubjectResult:
    config = job.config
    base_feature = config.base_features[0]
    data, rois = _load_feature(job)
    roi_lookup = np.asarray(job.atlas_info.roi_values)

    roi_rows: list[dict[str, Any]] = []
    payload: dict[Any, Any] = {}
    for stat_func, stat_name in zip(config.roi_stats, config.roi_stat_names):
        roi_values = calc_roi_statistics(data, rois, roi_lookup, stat_func)
        for roi, value in zip(job.atlas_info.node_labels, roi_values):
            roi_rows.append(
                {
                    "subject_id": job.subject_id,
                    "base_feature": base_feature,
                    "stat_name": stat_name,
                    "roi": roi,
                    "value": float(value),
                }
            )
        if config.return_results:
            key = job.subject_id if len(config.roi_stat_names) == 1 else (stat_name, job.subject_id)
            payload[key] = np.asarray(roi_values, dtype=float)

    return SubjectResult([], [], roi_rows, payload)


def _process_multiedge(job: SubjectJob) -> SubjectResult:
    config = job.config
    feature_cache = {
        base_feature: _load_feature(job, base_feature=base_feature)
        for base_feature in config.base_features
    }

    raw_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    payload: dict[Any, Any] = {}

    for weight_method in config.weight_methods:
        rows_for_weight: list[dict[str, Any]] = []
        for base_feature in config.base_features:
            data, rois = feature_cache[base_feature]
            graph = hiwenet.extract(
                data,
                rois,
                weight_method=weight_method,
                num_bins=config.num_bins,
                edge_range=config.edge_range_dict[base_feature],
                return_networkx_graph=True,
            )
            rows = _rows_from_graph(
                graph,
                job.subject_id,
                base_feature,
                weight_method,
                job.atlas_info.node_labels,
            )
            vector = _vector_from_rows(rows, job.atlas_info.node_labels)
            warn_nan(vector)
            rows_for_weight.extend(rows)
            raw_rows.extend(rows)
            if config.return_results:
                payload[(weight_method, base_feature, job.subject_id)] = vector

        grouped_weights: dict[tuple[str, str], list[float]] = defaultdict(list)
        for row in rows_for_weight:
            grouped_weights[(row["u"], row["v"])].append(row["weight"])

        for stat_func, stat_name in zip(config.summary_stats, config.summary_stat_names):
            for (u, v), weights in grouped_weights.items():
                summary_rows.append(
                    {
                        "subject_id": job.subject_id,
                        "weight_method": weight_method,
                        "summary_stat": stat_name,
                        "u": u,
                        "v": v,
                        "weight": float(stat_func(np.asarray(weights, dtype=float))),
                    }
                )

    return SubjectResult(raw_rows, summary_rows, [], payload)


def _process_job(job: SubjectJob) -> SubjectResult:
    if job.config.mode == "edges":
        return _process_edges(job)
    if job.config.mode == "roi-stats":
        return _process_roi_stats(job)
    if job.config.mode == "multiedge":
        return _process_multiedge(job)
    raise ValueError(f"Unknown graynet run mode: {job.config.mode}")


def _execute_jobs(jobs: tuple[SubjectJob, ...], num_procs: int) -> list[SubjectResult]:
    if num_procs <= 1 or len(jobs) <= 1:
        return [_process_job(job) for job in jobs]

    with ProcessPoolExecutor(max_workers=num_procs) as pool:
        return list(pool.map(_process_job, jobs))


def _write_metadata(run_dir: Path, config: RunConfig, atlas_info: AtlasInfo) -> None:
    metadata = {
        "mode": config.mode,
        "atlas_name": atlas_info.atlas_name,
        "atlas_spec": atlas_identifier(atlas_info.atlas_spec, atlas_info.atlas_name),
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


def _write_results(config: RunConfig, atlas_info: AtlasInfo, jobs: tuple[SubjectJob, ...]) -> tuple[dict[Any, Any], Path | None]:
    run_dir = build_run_dir(config)
    if run_dir is not None:
        _write_metadata(run_dir, config, atlas_info)

    raw_writer = ParquetBatchWriter(run_dir / RAW_EDGE_FILE_NAME, RAW_EDGE_SCHEMA) if run_dir else None
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
        for result in _execute_jobs(jobs, config.num_procs):
            if raw_writer:
                raw_writer.write_rows(result.raw_edge_rows)
            if summary_writer:
                summary_writer.write_rows(result.summary_edge_rows)
            if roi_writer:
                roi_writer.write_rows(result.roi_stat_rows)
            payload.update(result.return_payload)
    finally:
        for writer in writers:
            writer.close()

    return payload, run_dir


def _normalize_common(config: RunConfig) -> tuple[RunConfig, tuple[str, ...]]:
    subject_ids, _, _, _ = check_subjects(config.subject_ids)
    return (
        replace(
            config,
            input_dir=Path(config.input_dir),
            out_dir=Path(config.out_dir) if config.out_dir is not None else None,
            subject_ids=tuple(str(subject) for subject in subject_ids),
            num_procs=check_num_procs(config.num_procs),
        ),
        tuple(str(subject) for subject in subject_ids),
    )


def _resolve_run(config: RunConfig) -> tuple[RunConfig, AtlasInfo]:
    config, subject_ids = _normalize_common(config)

    if config.mode == "edges":
        base_feature = check_features(config.base_features[0])[0]
        check_params_single_edge(
            base_feature,
            config.input_dir,
            config.atlas,
            config.smoothing_param,
            config.node_size,
            config.out_dir,
            config.return_results,
        )
        atlas_info = resolve_atlas(base_feature, config.atlas, config.node_size)
        weight_methods, _, _, _ = check_weights(config.weight_methods)
        num_bins, edge_range = check_weight_params(config.num_bins, config.edge_range)
        return (
            replace(
                config,
                subject_ids=subject_ids,
                base_features=(base_feature,),
                weight_methods=tuple(weight_methods),
                atlas_name=atlas_info.atlas_name,
                num_bins=num_bins,
                edge_range=edge_range,
                num_procs=config.num_procs,
            ),
            atlas_info,
        )

    if config.mode == "roi-stats":
        base_feature = check_features(config.base_features[0])[0]
        check_params_single_edge(
            base_feature,
            config.input_dir,
            config.atlas,
            config.smoothing_param,
            config.node_size,
            config.out_dir,
            config.return_results,
        )
        atlas_info = resolve_atlas(base_feature, config.atlas, config.node_size)
        stat_funcs, stat_names, _, _, _ = check_stat_methods(config.roi_stats)
        return (
            replace(
                config,
                subject_ids=subject_ids,
                base_features=(base_feature,),
                atlas_name=atlas_info.atlas_name,
                roi_stats=tuple(stat_funcs),
                roi_stat_names=tuple(stat_names),
                num_procs=config.num_procs,
            ),
            atlas_info,
        )

    if config.mode == "multiedge":
        for feature in config.base_features:
            if feature in cfg.features_volumetric:
                raise NotImplementedError(
                    "Multi-edge networks are not yet supported for volumetric features."
                )
        check_params_multiedge(
            config.base_features,
            config.input_dir,
            config.atlas,
            config.smoothing_param,
            config.node_size,
            config.out_dir,
            config.return_results,
        )
        base_features = tuple(check_features(config.base_features))
        atlas_info = resolve_atlas(base_features[0], config.atlas, config.node_size)
        weight_methods, _, _, _ = check_weights(config.weight_methods)
        summary_stats, summary_names, _, _, _ = check_stat_methods(config.summary_stats)
        edge_range_dict = check_edge_range_dict(
            dict(config.edge_range_dict) if config.edge_range_dict is not None else None,
            base_features,
        )
        return (
            replace(
                config,
                subject_ids=subject_ids,
                base_features=base_features,
                weight_methods=tuple(weight_methods),
                atlas_name=atlas_info.atlas_name,
                num_bins=check_num_bins(config.num_bins),
                edge_range_dict=edge_range_dict,
                summary_stats=tuple(summary_stats),
                summary_stat_names=tuple(summary_names),
                num_procs=config.num_procs,
            ),
            atlas_info,
        )

    raise ValueError(f"Unknown graynet run mode: {config.mode}")


def _build_jobs(config: RunConfig, atlas_info: AtlasInfo) -> tuple[SubjectJob, ...]:
    return tuple(
        SubjectJob(subject_id=subject_id, config=config, atlas_info=atlas_info)
        for subject_id in config.subject_ids
    )


def run(config: RunConfig):
    resolved_config, atlas_info = _resolve_run(config)
    jobs = _build_jobs(resolved_config, atlas_info)
    payload, run_dir = _write_results(resolved_config, atlas_info, jobs)
    return payload if resolved_config.return_results else run_dir


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
    return run(
        RunConfig(
            mode="edges",
            input_dir=input_dir,
            out_dir=out_dir,
            subject_ids=subject_id_list,
            atlas=atlas,
            smoothing_param=smoothing_param,
            node_size=node_size,
            base_features=(base_feature,),
            weight_methods=tuple(weight_method_list),
            num_bins=num_bins,
            edge_range=edge_range,
            return_results=return_results,
            num_procs=num_procs,
        )
    )


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
    return run(
        RunConfig(
            mode="roi-stats",
            input_dir=input_dir,
            out_dir=out_dir,
            subject_ids=subject_id_list,
            atlas=atlas,
            smoothing_param=smoothing_param,
            node_size=node_size,
            base_features=(base_feature,),
            roi_stats=tuple(chosen_roi_stats),
            return_results=return_results,
            num_procs=num_procs,
        )
    )


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
    return run(
        RunConfig(
            mode="multiedge",
            input_dir=input_dir,
            out_dir=out_dir,
            subject_ids=subject_id_list,
            atlas=atlas,
            smoothing_param=smoothing_param,
            node_size=node_size,
            base_features=tuple(base_feature_list),
            weight_methods=tuple(weight_method_list),
            num_bins=num_bins,
            edge_range_dict=edge_range_dict,
            summary_stats=tuple(summary_stats),
            return_results=return_results,
            num_procs=num_procs,
            overwrite_results=overwrite_results,
        )
    )
