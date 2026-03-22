__all__ = ["extract_multiedge", "summarize_multigraph"]

import networkx as nx
import numpy as np

from graynet import config_graynet as cfg
from graynet.pipeline import run_multiedge


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
    return run_multiedge(
        subject_id_list,
        input_dir,
        base_feature_list=base_feature_list,
        weight_method_list=weight_method_list,
        summary_stats=summary_stats,
        num_bins=num_bins,
        edge_range_dict=edge_range_dict,
        atlas=atlas,
        smoothing_param=smoothing_param,
        node_size=node_size,
        out_dir=out_dir,
        return_results=return_results,
        overwrite_results=overwrite_results,
        num_procs=num_procs,
    )


def summarize_multigraph(multigraph, func_summary):
    summary_multigraph = nx.Graph()
    for u, v in multigraph.edges():
        all_weights = np.array(
            [edge_item["weight"] for _, edge_item in multigraph[u][v].items()],
            dtype=float,
        )
        summary_multigraph.add_edge(u, v, weight=float(func_summary(all_weights)))
    return summary_multigraph
