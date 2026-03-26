from graynet import config_graynet as cfg
from graynet.domain import GraynetJobConfig
from graynet.pipeline import run_graynet

__all__ = ["extract", "extract_multiedge", "roiwise_stats_indiv"]


def extract(
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


def extract_multiedge(
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


def roiwise_stats_indiv(
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
