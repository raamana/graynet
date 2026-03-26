from __future__ import annotations

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

import hiwenet
import networkx as nx
import numpy as np

from graynet import config_graynet as cfg
from graynet.atlas import atlas_identifier, resolve_atlas
from graynet.domain import AtlasInfo, GraynetJob, GraynetJobConfig, RunConfig, SubjectBatchResult
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


def _load_masked_feature(job: GraynetJob, base_feature: str | None = None):
    chosen_feature = base_feature or job.base_feature
    features = import_features(
        job.input_dir,
        [job.subject_id],
        chosen_feature,
        fwhm=job.smoothing_param,
        atlas=job.atlas_info.atlas_spec,
    )
    return mask_background_roi(
        features[job.subject_id],
        job.atlas_info.roi_labels,
        job.atlas_info.ignore_label,
    )


def _run_single_subject_edges(job: GraynetJob) -> SubjectBatchResult:
    raw_rows: list[dict[str, Any]] = []
    payload: dict[Any, Any] = {}
    data, rois = _load_masked_feature(job)

    for weight_method in job.weight_methods:
        graph = hiwenet.extract(
            data,
            rois,
            weight_method=weight_method,
            num_bins=job.num_bins,
            edge_range=job.edge_range,
            return_networkx_graph=True,
        )
        rows = _rows_from_graph(
            graph, job.subject_id, job.base_feature, weight_method, job.atlas_info.node_labels
        )
        vector = _vector_from_rows(rows, job.atlas_info.node_labels)
        warn_nan(vector)
        raw_rows.extend(rows)
        if job.return_results:
            payload[(weight_method, job.subject_id)] = vector

    return SubjectBatchResult(raw_rows, [], [], payload)


def _run_single_subject_roi_stats(job: GraynetJob) -> SubjectBatchResult:
    data, rois = _load_masked_feature(job)

    roi_stat_rows: list[dict[str, Any]] = []
    payload: dict[Any, Any] = {}
    unique_rois = job.atlas_info.roi_values
    roi_lookup = np.asarray(unique_rois)

    for stat_func, stat_name in zip(job.stat_funcs, job.stat_names):
        roi_values = calc_roi_statistics(data, rois, roi_lookup, stat_func)
        for roi, value in zip(job.atlas_info.node_labels, roi_values):
            roi_stat_rows.append(
                {
                    "subject_id": job.subject_id,
                    "base_feature": job.base_feature,
                    "stat_name": stat_name,
                    "roi": roi,
                    "value": float(value),
                }
            )
        if job.return_results:
            key = job.subject_id if len(job.stat_names) == 1 else (stat_name, job.subject_id)
            payload[key] = np.asarray(roi_values, dtype=float)

    return SubjectBatchResult([], [], roi_stat_rows, payload)


def _run_single_subject_multiedge(job: GraynetJob) -> SubjectBatchResult:
    feature_cache = {}
    for base_feature in job.base_features:
        feature_cache[base_feature] = _load_masked_feature(job, base_feature=base_feature)

    raw_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    payload: dict[Any, Any] = {}

    for weight_method in job.weight_methods:
        rows_for_weight: list[dict[str, Any]] = []
        for base_feature in job.base_features:
            data, rois = feature_cache[base_feature]
            graph = hiwenet.extract(
                data,
                rois,
                weight_method=weight_method,
                num_bins=job.num_bins,
                edge_range=job.edge_range_dict[base_feature],
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
            if job.return_results:
                payload[(weight_method, base_feature, job.subject_id)] = vector

        grouped_weights: dict[tuple[str, str], list[float]] = defaultdict(list)
        for row in rows_for_weight:
            grouped_weights[(row["u"], row["v"])].append(row["weight"])

        for stat_func, stat_name in zip(job.summary_stats, job.summary_stat_names):
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

    return SubjectBatchResult(raw_rows, summary_rows, [], payload)


def _process_job(job: GraynetJob) -> SubjectBatchResult:
    if job.mode == "edges":
        return _run_single_subject_edges(job)
    if job.mode == "roi-stats":
        return _run_single_subject_roi_stats(job)
    if job.mode == "multiedge":
        return _run_single_subject_multiedge(job)
    raise ValueError(f"Unknown Graynet job mode: {job.mode}")


def _collect_results(jobs: tuple[GraynetJob, ...], num_procs: int) -> list[SubjectBatchResult]:
    if num_procs <= 1 or len(jobs) <= 1:
        return [_process_job(job) for job in jobs]

    with ProcessPoolExecutor(max_workers=num_procs) as pool:
        return list(pool.map(_process_job, jobs))


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
    jobs: tuple[GraynetJob, ...],
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
        for batch in _collect_results(jobs, config.num_procs):
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


def _make_run_config(
    mode: str,
    input_dir,
    out_dir,
    subject_ids,
    atlas_info: AtlasInfo,
    smoothing_param,
    node_size,
    return_results,
    num_procs,
    *,
    base_features=(),
    weight_methods=(),
    num_bins=None,
    edge_range=None,
    edge_range_dict=None,
    summary_stats=(),
    summary_stat_names=(),
    roi_stats=(),
    roi_stat_names=(),
) -> RunConfig:
    return RunConfig(
        mode=mode,
        input_dir=Path(input_dir),
        out_dir=Path(out_dir) if out_dir is not None else None,
        subject_ids=tuple(str(subject) for subject in subject_ids),
        base_features=tuple(base_features),
        weight_methods=tuple(weight_methods),
        atlas_spec=atlas_info.atlas_spec,
        atlas_name=atlas_info.atlas_name,
        smoothing_param=smoothing_param,
        node_size=node_size,
        num_bins=num_bins,
        edge_range=edge_range,
        edge_range_dict=edge_range_dict,
        summary_stats=tuple(summary_stats),
        summary_stat_names=tuple(summary_stat_names),
        roi_stats=tuple(roi_stats),
        roi_stat_names=tuple(roi_stat_names),
        return_results=return_results,
        num_procs=num_procs,
    )


def _finalize_run(config: GraynetJobConfig) -> tuple[RunConfig, AtlasInfo, tuple[GraynetJob, ...]]:
    if config.mode == "edges":
        feature_name = check_features(config.base_features[0])[0]
        check_params_single_edge(
            feature_name,
            config.input_dir,
            config.atlas,
            config.smoothing_param,
            config.node_size,
            config.out_dir,
            config.return_results,
        )
        atlas_info = resolve_atlas(feature_name, config.atlas, config.node_size)
        subject_ids, _, _, _ = check_subjects(config.subject_ids)
        weight_methods, _, _, _ = check_weights(config.weight_methods)
        num_bins, edge_range = check_weight_params(config.num_bins, config.edge_range)
        num_procs = check_num_procs(config.num_procs)
        config = _make_run_config(
            "edges",
            config.input_dir,
            config.out_dir,
            subject_ids,
            atlas_info,
            config.smoothing_param,
            config.node_size,
            config.return_results,
            num_procs,
            base_features=(feature_name,),
            weight_methods=weight_methods,
            num_bins=num_bins,
            edge_range=edge_range,
        )
        jobs = tuple(
            GraynetJob(
                mode="edges",
                subject_id=subject_id,
                input_dir=config.input_dir,
                atlas_info=atlas_info,
                base_feature=feature_name,
                weight_methods=config.weight_methods,
                smoothing_param=config.smoothing_param,
                num_bins=num_bins,
                edge_range=edge_range,
                return_results=config.return_results,
            )
            for subject_id in config.subject_ids
        )
    elif config.mode == "roi-stats":
        feature_name = check_features(config.base_features[0])[0]
        check_params_single_edge(
            feature_name,
            config.input_dir,
            config.atlas,
            config.smoothing_param,
            config.node_size,
            config.out_dir,
            config.return_results,
        )
        stat_funcs, stat_names, _, _, _ = check_stat_methods(config.roi_stats)
        atlas_info = resolve_atlas(feature_name, config.atlas, config.node_size)
        subject_ids, _, _, _ = check_subjects(config.subject_ids)
        num_procs = check_num_procs(config.num_procs)
        config = _make_run_config(
            "roi-stats",
            config.input_dir,
            config.out_dir,
            subject_ids,
            atlas_info,
            config.smoothing_param,
            config.node_size,
            config.return_results,
            num_procs,
            base_features=(feature_name,),
            roi_stats=stat_funcs,
            roi_stat_names=stat_names,
        )
        jobs = tuple(
            GraynetJob(
                mode="roi-stats",
                subject_id=subject_id,
                input_dir=config.input_dir,
                atlas_info=atlas_info,
                base_feature=feature_name,
                smoothing_param=config.smoothing_param,
                stat_funcs=config.roi_stats,
                stat_names=config.roi_stat_names,
                return_results=config.return_results,
            )
            for subject_id in config.subject_ids
        )
    elif config.mode == "multiedge":
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
        features = tuple(check_features(config.base_features))
        subject_ids, _, _, _ = check_subjects(config.subject_ids)
        weight_methods, _, _, _ = check_weights(config.weight_methods)
        stat_funcs, stat_names, _, _, _ = check_stat_methods(config.summary_stats)
        num_bins = check_num_bins(config.num_bins)
        given_ranges = dict(config.edge_range_dict) if config.edge_range_dict is not None else None
        edge_range_dict = check_edge_range_dict(given_ranges, features)
        num_procs = check_num_procs(config.num_procs)
        atlas_info = resolve_atlas(features[0], config.atlas, config.node_size)
        config = _make_run_config(
            "multiedge",
            config.input_dir,
            config.out_dir,
            subject_ids,
            atlas_info,
            config.smoothing_param,
            config.node_size,
            config.return_results,
            num_procs,
            base_features=features,
            weight_methods=weight_methods,
            num_bins=num_bins,
            edge_range_dict=edge_range_dict,
            summary_stats=stat_funcs,
            summary_stat_names=stat_names,
        )
        jobs = tuple(
            GraynetJob(
                mode="multiedge",
                subject_id=subject_id,
                input_dir=config.input_dir,
                atlas_info=atlas_info,
                base_features=config.base_features,
                weight_methods=config.weight_methods,
                smoothing_param=config.smoothing_param,
                num_bins=num_bins,
                edge_range_dict=edge_range_dict,
                summary_stats=config.summary_stats,
                summary_stat_names=config.summary_stat_names,
                return_results=config.return_results,
            )
            for subject_id in config.subject_ids
        )
    else:
        raise ValueError(f"Unknown graynet run mode: {config.mode}")

    return config, atlas_info, jobs


def run_graynet(config: GraynetJobConfig):
    finalized_config, atlas_info, jobs = _finalize_run(config)
    payload, run_dir = _run_with_writers(finalized_config, atlas_info, jobs)
    return payload if finalized_config.return_results else run_dir


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
    return run_graynet(
        GraynetJobConfig(
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
    return run_graynet(
        GraynetJobConfig(
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
    return run_graynet(
        GraynetJobConfig(
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
