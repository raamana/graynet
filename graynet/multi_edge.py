__all__ = ['extract', ]

import collections
import os
import sys
import warnings
import traceback
import logging
from os.path import join as pjoin, exists as pexists
from multiprocessing import Manager, Pool
from functools import partial

import numpy as np
import networkx as nx
import hiwenet

from sys import version_info

if version_info.major > 2:
    from graynet import parcellate
    from graynet import freesurfer
    from graynet import config_graynet as cfg
    from graynet import run_workflow as single_edge
else:
    raise NotImplementedError('graynet supports only Python 2.7 or 3+. Upgrade to Python 3+ is recommended.')

def extract(subject_id_list,
            input_dir,
            base_feature_list=cfg.default_feature_list_multi_edge,
            weight_method_list=cfg.default_weight_method,
            num_bins=cfg.default_num_bins,
            edge_range=cfg.default_edge_range,
            atlas=cfg.default_atlas,
            smoothing_param=cfg.default_smoothing_param,
            node_size=cfg.default_node_size,
            out_dir=None,
            return_results=False,
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

    num_bins : int
        Number of histogram bins to use when computing pair-wise weights based on histogram distance. Default : 25

    edge_range : tuple or list
        The range of edges (two finite values) within which to build the histogram e.g. ``--edge_range 0 5``.
        This can be helpful (and important) to ensure correspondence across multiple invocations of graynet (e.g. for different subjects), in terms of range across all bins as well as individual bin edges.

        Default :

            - ( 0.0, 5.0) for ``freesurfer_thickness`` and
            - (-0.3, 0.3) for ``freesurfer_curv``.

    atlas : str
        Name of the atlas whose parcellation to be used.
        Choices for cortical parcellation: ['FSAVERAGE', 'GLASSER2016'], which are primary cortical.
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
    check_params_multiedge(base_feature_list, input_dir, atlas, smoothing_param, node_size, out_dir, return_results)
    atlas = single_edge.check_atlas(atlas)

    subject_id_list, num_subjects, max_id_width, nd_id = single_edge.check_subjects(subject_id_list)

    num_bins, edge_range = single_edge.check_weight_params(num_bins, edge_range)
    weight_method_list, num_weights, max_wtname_width, nd_wm = single_edge.check_weights(weight_method_list)

    num_procs = single_edge.check_num_procs(num_procs)
    pretty_print_options = (max_id_width, nd_id, num_weights, max_wtname_width, nd_wm)

    # roi_labels, ctx_annot = parcellate.freesurfer_roi_labels(atlas)
    # uniq_rois, roi_size, num_nodes = roi_info(roi_labels)
    uniq_rois, centroids, roi_labels = parcellate.roi_labels_centroids(atlas)

    print('\nProcessing {} features resampled to {} atlas,'
          ' smoothed at {} with node size {}'.format(base_feature_list, atlas,
                                                     smoothing_param, node_size))

    if not return_results:
        if out_dir is None:
            raise ValueError('When return_results=False, out_dir must be specified to be able to save the results.')
        if not pexists(out_dir):
            os.mkdir(out_dir)

    chunk_size = int(np.ceil(num_subjects / num_procs))
    with Manager():
        partial_func_extract = partial(per_subject_multi_edge, input_dir, base_feature_list,
                                       roi_labels, centroids,
                                       weight_method_list,
                                       atlas, smoothing_param, node_size,
                                       num_bins, edge_range,
                                       out_dir, return_results, pretty_print_options)
        with Pool(processes=num_procs) as pool:
            edge_weights_list_dicts = pool.map(partial_func_extract, subject_id_list, chunk_size)

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
                           weight_method_list,
                           atlas, smoothing_param, node_size,
                           num_bins, edge_range,
                           out_dir, return_results, pretty_print_options,
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

    func_summary = np.nanmedian

    for ww, weight_method in enumerate(weight_method_list):

        multigraph = nx.MultiGraph(weight_method=weight_method,
                                   subject_id=subject)

        for base_feature in base_feature_list:
            try:
                features = single_edge.import_features(input_dir,
                                                       [subject, ],
                                                       base_feature,
                                                       fwhm=smoothing_param,
                                                       atlas=atlas)
                # features[base_feature] is another dict keyed in by subject id
            except:
                traceback.print_exc()
                warnings.warn('Unable to read {} features'
                              ' for {}\n Skipping it.'.format(base_feature, subject), UserWarning)
                return

            data, rois = single_edge.mask_background_roi(features[subject], roi_labels,
                                                         parcellate.null_roi_name)

            # unique stamp for each subject and weight
            expt_id = single_edge.stamp_expt_weight(base_feature, atlas, smoothing_param, node_size, weight_method)
            sys.stdout.write('\nProcessing id {:{id_width}} -- weight {:{wtname_width}} ({:{nd_wm}}/{:{nd_wm}}) :'.format(subject, weight_method, ww + 1, num_weights, nd_id=nd_id, nd_wm=nd_wm, id_width=max_id_width, wtname_width=max_wtname_width))

            # actual computation of pair-wise features
            try:
                unigraph = hiwenet.extract(data,
                                           rois,
                                           weight_method=weight_method,
                                           num_bins=num_bins,
                                           edge_range=edge_range,
                                           return_networkx_graph=True)

                # retrieving edge weights
                weight_vec = np.array(list(nx.get_edge_attributes(unigraph, 'weight').values()))
                single_edge.warn_nan(weight_vec)
                if return_results:
                    edge_weights_all[(weight_method, base_feature, subject)] = weight_vec

            except (RuntimeError, RuntimeWarning) as runexc:
                print(runexc)
            except KeyboardInterrupt:
                print('Exiting on keyborad interrupt! \n'
                      'Abandoning the remaining processing ')
                sys.exit(1)
            except:
                print('Unable to extract {} weights for {} for {}'.format(weight_method, base_feature, subject))
                traceback.print_exc()

            print('Done.')

            # adding edges/weights from each feature to a multigraph
            # this also encodes the sources
            for u, v in unigraph.edges():
                multigraph.add_edge(u, v,
                                    weight=unigraph[u][v]['weight'],
                                    base_feature=base_feature)

        # adding position info to nodes (for visualization later)
        add_nodal_positions(multigraph, centroids)

        # creating single graph with a summary edge weight (like median)
        summary_multigraph = summarize_multigraph(multigraph, func_summary)
        add_nodal_positions(summary_multigraph, centroids)

        # saving to disk
        try:
            save_summary_graph(summary_multigraph, out_dir, subject,
                               expt_id, summary_descr=func_summary.__name__)
        except:
            raise IOError('Unable to save the graph to:\n{}'.format(out_dir))

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


def check_params_multiedge(base_feature_list, input_dir, atlas, smoothing_param,
                           node_size, out_dir, return_results):
    """Validation of parameters and appropriate type casting if necessary."""

    given = set(base_feature_list)
    allowed = set(cfg.base_feature_list)
    if not given.issubset(allowed):
        unrecog_methods = given.difference(allowed)
        raise NotImplementedError('Methods unrecognized: \n {}'.format(unrecog_methods))

    if atlas.upper() not in parcellate.atlas_list:
        raise ValueError('Invalid atlas choice. Use one of {}'.format(parcellate.atlas_list))

    if not pexists(input_dir):
        raise IOError('Input directory at {} does not exist.'.format(input_dir))

    if out_dir is None and return_results is False:
        raise ValueError('Results are neither saved to disk or being received when returned.\n'
                         'Specify out_dir (not None) or make return_results=True')

    if out_dir is not None and not pexists(out_dir):
        os.mkdir(out_dir)

    # no checks on subdivison size yet, as its not implemented

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
            out_file_name = '{}_{}_multigraph_graynet.graphml'.format(str_suffix,summary_descr)
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