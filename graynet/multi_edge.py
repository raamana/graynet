import graynet.utils

__all__ = ['extract_multiedge', 'summarize_multigraph']

import os
import sys
import warnings
import traceback
from os.path import join as pjoin, exists as pexists, isfile, realpath, getsize
from multiprocessing import Manager, Pool
from functools import partial

import numpy as np
import networkx as nx
import hiwenet

from sys import version_info

if version_info.major > 2:
    from graynet.utils import stamp_expt_multiedge, check_params_multiedge, make_output_path_graph, \
        save_graph, check_subjects, check_stat_methods, check_num_bins, check_weights, \
        check_num_procs, check_atlas, check_edge_range_dict, mask_background_roi, warn_nan, \
        stamp_expt_weight, import_features, save_per_subject_graph
    from graynet import parcellate
    from graynet import config_graynet as cfg
    from graynet import run_workflow as single_edge
else:
    raise NotImplementedError(
        'graynet supports only Python 2.7 or 3+. Upgrade to Python 3+ is recommended.')


def extract_multiedge(subject_id_list,
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
                      num_procs=cfg.default_num_procs):
    """
    Extracts weighted networks (matrix of pair-wise ROI distances) based on multiple gray matter features based on Freesurfer processing.

    Parameters
    ----------
    subject_id_list : str or list
         must be path to a file containing subject IDs, or a list of subject IDs
    input_dir : str
        Path to the input directory where features can be read.
        For example, this can be Freesurfer's SUBJECTS_DIR, where output processing is stored.
        Or another directory with a structure that graynet can parse.
    base_feature_list : list
        Set of features that drive the different edges between the pair of ROIs.

        For example, if you choose thickness and pial_curv, each pair of ROIs will have two edges.

        This multi-edge network can be turned into a single network based on averaging weights from different individual networks.

    weight_method : string(s), optional
        Type of distance (or metric) to compute between the pair of histograms.

        It must be one of the following methods:

        - 'chebyshev'
        - 'chebyshev_neg'
        - 'chi_square'
        - 'correlate'
        - 'correlate_1'
        - 'cosine'
        - 'cosine_1'
        - 'cosine_2'
        - 'cosine_alt'
        - 'euclidean'
        - 'fidelity_based'
        - 'histogram_intersection'
        - 'histogram_intersection_1'
        - 'jensen_shannon'
        - 'kullback_leibler'
        - 'manhattan'
        - 'minowski'
        - 'noelle_1'
        - 'noelle_2'
        - 'noelle_3'
        - 'noelle_4'
        - 'noelle_5'
        - 'relative_bin_deviation'
        - 'relative_deviation'

        Note only the following are *metrics*:

        - 'manhattan'
        - 'minowski'
        - 'euclidean'
        - 'noelle_2'
        - 'noelle_4'
        - 'noelle_5'

        The following are *semi- or quasi-metrics*:

        - 'kullback_leibler'
        - 'jensen_shannon'
        - 'chi_square'
        - 'chebyshev'
        - 'cosine_1'
        - 'chebyshev_neg'
        - 'correlate_1'
        - 'histogram_intersection_1'
        - 'relative_deviation'
        - 'relative_bin_deviation'
        - 'noelle_1'
        - 'noelle_3'

        The following are  classified to be similarity functions:

        - 'histogram_intersection'
        - 'correlate'
        - 'cosine'
        - 'cosine_2'
        - 'cosine_alt'
        - 'fidelity_based'

        *Default* choice: 'manhattan'.

    summary_stats : list of str
        A string, or list of strings, each representing a method (like 'median', 'prod' or 'max'),
        to compute a summay statistic from the array of multiple weights computed.

        This must be available as a member of numpy or scipy.stats.

    num_bins : int
        Number of histogram bins to use when computing pair-wise weights based on histogram distance. Default : 25

    edge_range_dict : tuple or list
        The range of edges (two finite values) within which to build the histogram e.g. ``--edge_range 0 5``.
        This can be helpful (and important) to ensure correspondence across multiple invocations of graynet (e.g. for different subjects), in terms of range across all bins as well as individual bin edges.

        Default :

            - ( 0.0, 5.0) for ``freesurfer_thickness`` and
            - (-0.3, 0.3) for ``freesurfer_curv``.

    atlas : str
        Name of the atlas whose parcellation to be used.
        Choices for cortical parcellation: ['fsaverage', 'glasser2016'], which are primary cortical.
        Volumetric whole-brain atlases will be added soon.

    smoothing_param : scalar
        Smoothing parameter, which could be fwhm for Freesurfer cortical features,
        or another relevant for the chosen base_feature_list.
        Default: assumed as fwhm=10mm for the default feature choice 'thickness'

    node_size : scalar, optional
        Parameter to indicate the size of the ROIs, subparcels or patches, depending on type of atlas or feature.
        This feature is not implemented yet, just a placeholder and to enable default computation.

    out_dir : str, optional
        Path to output directory to store results.
        Default: None, results are returned, but not saved to disk.
        If this is None, return_results must be true.

    return_results : bool
        Flag to indicate whether to return the results to be returned.
        This flag helps to reduce the memory requirements, when the number of nodes in a parcellation or
        the number of subjects or weight methods are large, as it doesn't retain results for all combinations,
        when running from commmand line interface (or HPC). Default: False
        If this is False, out_dir must be specified to save the results to disk.

    overwrite_results : bool
        Flag to request overwriting of existing results, in case of reruns/failed jobs. By default, if the expected output file exists and is of non-zero size, its computation is skipped (assuming the file is complete, usable and not corrupted).

    num_procs : int
        Number of parallel processes to use to speed up computation.

    Returns
    -------
    edge_weights_all : dict, None
        If return_results is True, this will be a dictionary keyed in by a tuple: (weight method, subject_ID)
        The value of each edge_weights_all[(weight method, subject_ID)] is
        a numpy array of length p = k*(k-1)/2, with k = number of nodes in the atlas parcellation.
        If return_results is False, this will be None, which is the default.
    """

    # All the checks must happen here, as this is key function in the API
    check_params_multiedge(base_feature_list, input_dir, atlas, smoothing_param, node_size, out_dir,
                           return_results)
    atlas = check_atlas(atlas)

    subject_id_list, num_subjects, max_id_width, nd_id = check_subjects(subject_id_list)

    num_bins = check_num_bins(num_bins)
    edge_range_dict = check_edge_range_dict(edge_range_dict, base_feature_list)
    weight_method_list, num_weights, max_wtname_width, nd_wm = check_weights(
        weight_method_list)

    # validating the choice and getting a callable
    summary_stats, summary_stat_names, _, _, _ = check_stat_methods(summary_stats)

    num_procs = check_num_procs(num_procs)
    pretty_print_options = (max_id_width, nd_id, num_weights, max_wtname_width, nd_wm)

    # roi_labels, ctx_annot = parcellate.freesurfer_roi_labels(atlas)
    # uniq_rois, roi_size, num_nodes = roi_info(roi_labels)
    uniq_rois, centroids, roi_labels = parcellate.roi_labels_centroids(atlas)

    print('\nProcessing {} features resampled to {} atlas,'
          ' smoothed at {} with node size {}'.format(base_feature_list, atlas,
                                                     smoothing_param, node_size))

    if not return_results:
        if out_dir is None:
            raise ValueError(
                'When return_results=False, out_dir must be specified to be able to save the results.')
        if not pexists(out_dir):
            os.mkdir(out_dir)

    partial_func_extract = partial(per_subject_multi_edge, input_dir, base_feature_list,
                                   roi_labels, centroids,
                                   weight_method_list, summary_stats, summary_stat_names,
                                   atlas, smoothing_param, node_size,
                                   num_bins, edge_range_dict,
                                   out_dir, return_results, overwrite_results, pretty_print_options)
    if num_procs > 1:
        chunk_size = int(np.ceil(num_subjects / num_procs))
        with Manager():
            with Pool(processes=num_procs) as pool:
                edge_weights_list_dicts = pool.map(partial_func_extract, subject_id_list,
                                                   chunk_size)
    else:
        # reverting to sequential processing
        edge_weights_list_dicts = [partial_func_extract(subject=sub_id) for sub_id in
                                   subject_id_list]

    if return_results:
        edge_weights_all = dict()
        for combo in edge_weights_list_dicts:
            # each element from output of parallel loop is a dict keyed in by {subject, weight)
            edge_weights_all.update(combo)
    else:
        edge_weights_all = None

    print('\ngraynet computation done.')
    return edge_weights_all


def per_subject_multi_edge(input_dir, base_feature_list, roi_labels, centroids,
                           weight_method_list, summary_stats, summary_stat_names,
                           atlas, smoothing_param, node_size,
                           num_bins, edge_range_dict,
                           out_dir, return_results, overwrite_results, pretty_print_options,
                           subject=None):  # purposefully leaving it last to enable partial function creation
    """
    Extracts give set of weights for one subject.
    """

    if subject is None:
        return

    if return_results:
        edge_weights_all = dict()
    else:
        edge_weights_all = None

    max_id_width, nd_id, num_weights, max_wtname_width, nd_wm = pretty_print_options

    for ww, weight_method in enumerate(weight_method_list):

        expt_id_multi = stamp_expt_multiedge(base_feature_list, atlas, smoothing_param, node_size,
                                             weight_method)
        out_path_multigraph = make_output_path_graph(out_dir, subject,
                                                     [expt_id_multi, 'multigraph'])
        # skipping the computation if the file exists already
        if not overwrite_results and isfile(out_path_multigraph) and getsize(
                out_path_multigraph) > 0:
            print('\nMultigraph exists already at\n\t{}\n'
                  ' skipping its computation!'.format(out_path_multigraph))
            multigraph = None  # signal to re-read
        else:
            multigraph = nx.MultiGraph()

            for base_feature in base_feature_list:
                # # TODO refactor
                # unigraph, weight_vec = compute_unigraph(input_dir, subject, base_feature, weight_method, roi_labels,
                #                                         atlas, smoothing_param, node_size, centroids,
                #                                         num_bins, edge_range_dict,
                #                                         out_dir, overwrite_results, pretty_print_options)
                # if return_results:
                #     edge_weights_all[(weight_method, base_feature, subject)] = weight_vec

                try:
                    features = import_features(input_dir,
                                               [subject, ],
                                               base_feature,
                                               fwhm=smoothing_param,
                                               atlas=atlas)

                except:
                    traceback.print_exc()
                    warnings.warn('Unable to read {} features'
                                  ' for {}\n Skipping it.'.format(base_feature, subject),
                                  UserWarning)
                    return

                data, rois = mask_background_roi(features[subject], roi_labels,
                                                               cfg.null_roi_name)

                # unique stamp for each subject and weight
                expt_id_single = stamp_expt_weight(base_feature, atlas, smoothing_param,
                                                                 node_size, weight_method)
                sys.stdout.write('\nProcessing id {:{id_width}} --'
                                 ' weight {:{wtname_width}} ({:{nd_wm}}/{:{nd_wm}})'
                                 ' :'.format(subject, weight_method, ww + 1, num_weights,
                                             nd_id=nd_id, nd_wm=nd_wm, id_width=max_id_width,
                                             wtname_width=max_wtname_width))

                # actual computation of pair-wise features
                try:
                    unigraph = hiwenet.extract(data,
                                               rois,
                                               weight_method=weight_method,
                                               num_bins=num_bins,
                                               edge_range=edge_range_dict[base_feature],
                                               return_networkx_graph=True)

                    # retrieving edge weights
                    weight_vec = np.array(list(nx.get_edge_attributes(unigraph, 'weight').values()))
                    warn_nan(weight_vec)
                    if return_results:
                        edge_weights_all[(weight_method, base_feature, subject)] = weight_vec

                except (RuntimeError, RuntimeWarning) as runexc:
                    print(runexc)
                except KeyboardInterrupt:
                    print('Exiting on keyborad interrupt! \n'
                          'Abandoning the remaining processing ')
                    sys.exit(1)
                except:
                    print('Unable to extract {} weights for {} for {}'.format(weight_method,
                                                                              base_feature,
                                                                              subject))
                    traceback.print_exc()

                print('Done.')

                # TODO consider extracting some network features upon user request.

                add_nodal_positions(unigraph, centroids)
                save_per_subject_graph(unigraph, out_dir, subject, expt_id_single)

                # adding edges/weights from each feature to a multigraph
                # this also encodes the sources
                for u, v in unigraph.edges():
                    multigraph.add_edge(u, v,
                                        weight=unigraph[u][v]['weight'],
                                        base_feature=base_feature)

            # adding position info to nodes (for visualization later)
            add_nodal_positions(multigraph, centroids)
            save_graph(multigraph, out_path_multigraph, 'multi-edge')

        for stat_func, stat_name in zip(summary_stats, summary_stat_names):
            # creating single graph with a summary edge weight (like median)
            out_path_summary = make_output_path_graph(out_dir, subject,
                                                      [expt_id_multi, stat_name, 'multigraph'])
            if not overwrite_results and isfile(out_path_summary) and getsize(out_path_summary) > 0:
                print(
                    'Summary {} of multigraph exists already at\n\t{}\n skipping its computation!'.format(
                        stat_name, out_path_summary))
            else:
                if multigraph is None:
                    multigraph = nx.read_graphml(out_path_multigraph)

                try:
                    summary_multigraph = summarize_multigraph(multigraph, stat_func)
                    add_nodal_positions(summary_multigraph, centroids)
                    save_graph(summary_multigraph, out_path_summary, '{} summary'.format(stat_name))
                except:
                    print('Summary {} could not be computed - skipping!'.format(stat_name))
                    traceback.print_exc()

    return edge_weights_all


def summarize_multigraph(multigraph, func_summary):
    "Creating single graph with a summary edge weight (like median)"

    summary_multigraph = nx.Graph()
    for u, v in multigraph.edges():
        # looping through parallel edges and obtaining their weights.
        all_weights = np.array([edge_item['weight'] for idx, edge_item in multigraph[u][v].items()])
        summary_weight = float(func_summary(all_weights))  # float needed due to graphml limitation
        summary_multigraph.add_edge(u, v, weight=summary_weight)

    return summary_multigraph


def add_nodal_positions(graph, centroids):
    "Adds the x, y, z attributes to each node in graph."

    # adding position info to nodes (for visualization later)
    for roi in centroids:
        graph.node[roi]['x'] = float(centroids[roi][0])
        graph.node[roi]['y'] = float(centroids[roi][1])
        graph.node[roi]['z'] = float(centroids[roi][2])

    return


def save_summary_graph(graph, out_dir, subject,
                       str_suffix=None,
                       summary_descr='summary'):
    "Saves the features to disk."

    if out_dir is not None:
        # get outpath returned from hiwenet, based on dist name and all other parameters
        # choose out_dir name  based on dist name and all other parameters
        out_subject_dir = pjoin(out_dir, subject)
        if not pexists(out_subject_dir):
            os.mkdir(out_subject_dir)

        if str_suffix is not None:
            out_file_name = '{}_{}_multigraph_graynet.graphml'.format(str_suffix, summary_descr)
        else:
            out_file_name = '_{}_multigraph_graynet.graphml'.format(summary_descr)

        out_weights_path = pjoin(out_subject_dir, out_file_name)

        try:
            nx.info(graph)
            nx.write_graphml(graph, out_weights_path, encoding='utf-8')
            print('\nSaved the summary multi-graph to \n{}'.format(out_weights_path))
        except:
            print('\nUnable to save summary multi-graph to \n{}'.format(out_weights_path))
            traceback.print_exc()

    return
