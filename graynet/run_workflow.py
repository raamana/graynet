from graynet.utils import check_atlas, check_num_procs, check_stat_methods, warn_nan, check_subjects, check_weights, \
    check_weight_params, check_params_single_edge

__all__ = ['extract', 'roiwise_stats_indiv', 'roi_info', 'cli_run']

import os
import sys
import argparse
import warnings
import traceback
import logging
from os.path import join as pjoin, exists as pexists
from multiprocessing import Manager, Pool
from functools import partial

import hiwenet
import numpy as np
import networkx as nx

from sys import version_info

if version_info.major == 2 and version_info.minor == 7:
    import freesurfer
    import parcellate
    import config_graynet as cfg
elif version_info.major > 2:
    from graynet import parcellate
    from graynet import freesurfer
    from graynet import config_graynet as cfg
    from graynet import utils
else:
    raise NotImplementedError('graynet supports only Python 2.7 or 3+. Upgrade to Python 3+ is recommended.')

np.seterr(divide='ignore', invalid='ignore')

# logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

def extract(subject_id_list,
            input_dir,
            base_feature=cfg.default_feature_single_edge,
            weight_method_list=cfg.default_weight_method,
            num_bins=cfg.default_num_bins,
            edge_range=cfg.default_edge_range,
            atlas=cfg.default_atlas, smoothing_param=cfg.default_smoothing_param,
            node_size=cfg.default_node_size,
            out_dir=None, return_results=False, num_procs=cfg.default_num_procs):
    """
    Extracts weighted networks (matrix of pair-wise ROI distances) from gray matter features based on Freesurfer processing.

    Parameters
    ----------
    subject_id_list : str or list
         must be path to a file containing subject IDs, or a list of subject IDs
    input_dir : str
        Path to the input directory where features can be read.
        For example, this can be Freesurfer's SUBJECTS_DIR, where output processing is stored.
        Or another directory with a structure that graynet can parse.
    base_feature : str
        Specific type of feature to read for each subject from the input directory.

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
        or another relevant for the chosen base_feature.
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
    check_params_single_edge(base_feature, input_dir, atlas, smoothing_param, node_size, out_dir, return_results)
    atlas = check_atlas(atlas)

    subject_id_list, num_subjects, max_id_width, nd_id = check_subjects(subject_id_list)

    num_bins, edge_range = check_weight_params(num_bins, edge_range)
    weight_method_list, num_weights, max_wtname_width, nd_wm = check_weights(weight_method_list)

    num_procs = check_num_procs(num_procs)
    pretty_print_options = (max_id_width, nd_id, num_weights, max_wtname_width, nd_wm)

    # roi_labels, ctx_annot = parcellate.freesurfer_roi_labels(atlas)
    # uniq_rois, roi_size, num_nodes = roi_info(roi_labels)
    uniq_rois, centroids, roi_labels = parcellate.roi_labels_centroids(atlas)

    print('\nProcessing {} features resampled to {} atlas,'
          ' smoothed at {} with node size {}'.format(base_feature, atlas, smoothing_param, node_size))

    if not return_results:
        if out_dir is None:
            raise ValueError('When return_results=False, out_dir must be specified to be able to save the results.')
        if not pexists(out_dir):
            os.mkdir(out_dir)

    chunk_size = int(np.ceil(num_subjects/num_procs))
    with Manager():
        partial_func_extract = partial(extract_per_subject, input_dir, base_feature, roi_labels, centroids, weight_method_list,
                                       atlas, smoothing_param, node_size, num_bins, edge_range, out_dir,
                                       return_results, pretty_print_options)
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


def extract_per_subject(input_dir, base_feature, roi_labels, centroids,
                        weight_method_list,
                        atlas, smoothing_param, node_size,
                        num_bins, edge_range,
                        out_dir, return_results, pretty_print_options,
                        subject=None): # purposefully leaving it last to enable partial function creation
    """
    Extracts give set of weights for one subject.

    Parameters
    ----------
    subject
    input_dir
    base_feature
    roi_labels
    weight_method_list
    atlas
    smoothing_param
    node_size
    num_bins
    edge_range
    out_dir
    return_results
    pretty_print_options

    Returns
    -------

    """

    if subject is None:
        return

    try:
        features = import_features(input_dir,
                                   [subject, ],
                                   base_feature,
                                   fwhm=smoothing_param,
                                   atlas=atlas)
    except:
        traceback.print_exc()
        warnings.warn('Unable to read {} features'
                      ' for {}\n Skipping it.'.format(base_feature, subject), UserWarning)

    data, rois = mask_background_roi(features[subject], roi_labels, parcellate.null_roi_name)

    max_id_width, nd_id, num_weights, max_wtname_width, nd_wm = pretty_print_options

    if return_results:
        edge_weights_all = dict()
    else:
        edge_weights_all = None

    for ww, weight_method in enumerate(weight_method_list):
        # unique stamp for each subject and weight
        expt_id = stamp_expt_weight(base_feature, atlas, smoothing_param, node_size, weight_method)
        sys.stdout.write('\nProcessing id {:{id_width}} -- weight {:{wtname_width}} ({:{nd_wm}}/{:{nd_wm}})'
                         ' :'.format(subject, weight_method, ww + 1, num_weights, nd_id=nd_id, nd_wm=nd_wm,
                                     id_width=max_id_width, wtname_width=max_wtname_width))

        # actual computation of pair-wise features
        try:
            graph = hiwenet.extract(data,
                                       rois,
                                       weight_method=weight_method,
                                       num_bins=num_bins,
                                       edge_range=edge_range,
                                       return_networkx_graph=True)

            # retrieving edge weights
            weight_vec = np.array(list(nx.get_edge_attributes(graph, 'weight').values()))
            warn_nan(weight_vec)
            # weight_vec = get_triu_handle_inf_nan(edge_weights)

            # adding position info to nodes (for visualization later)
            for roi in centroids:
                graph.node[roi]['x'] = float(centroids[roi][0])
                graph.node[roi]['y'] = float(centroids[roi][1])
                graph.node[roi]['z'] = float(centroids[roi][2])

            if return_results:
                edge_weights_all[(weight_method, subject)] = weight_vec

            # saving to disk
            try:
                save(weight_vec, out_dir, subject, expt_id)
                save_graph(graph, out_dir, subject, expt_id)
            except:
                raise IOError('Unable to save the network/vectorized weights to:\n{}'.format(out_dir))

        except (RuntimeError, RuntimeWarning) as runexc:
            print(runexc)
        except KeyboardInterrupt:
            print('Exiting on keyborad interrupt! \n'
                  'Abandoning the remaining processing for {} weights:\n'
                  '{}.'.format(num_weights - ww, weight_method_list[ww:]))
            sys.exit(1)
        except:
            print('Unable to extract {} features for {}'.format(weight_method, subject))
            traceback.print_exc()

        sys.stdout.write('Done.')


    return  edge_weights_all


def roiwise_stats_indiv(subject_id_list, input_dir,
                        base_feature=cfg.default_feature_single_edge,
                        chosen_roi_stats=cfg.default_roi_statistic,
                        atlas=cfg.default_atlas, smoothing_param=cfg.default_smoothing_param,
                        node_size=cfg.default_node_size,
                        out_dir=None, return_results=False):
    """
    Computes the chosen summary statistics within each ROI.
    These summary stats (such as median) can serve as a baseline for network-level values produced by graynet.

    Options for summary statistics include 'median', 'entropy', 'kurtosis' and
    any other appropriate summary statistics listed under scipy.stats:
    https://docs.scipy.org/doc/scipy/reference/stats.html#statistical-functions

    Parameters
    ----------
    subject_id_list : str or list
        must be path to a file containing subject IDs, or a list of subject IDs

    input_dir : str
        Path to the input directory where features can be read.
        For example, this can be Freesurfer's SUBJECTS_DIR, where output processing is stored.
        Or another directory with a structure that graynet can parse.

    base_feature : str
        Specific type of feature to read for each subject from the input directory.

    chosen_roi_stats : list of str or callable
        If requested, graynet will compute chosen summary statistics (such as median) within each ROI of the chosen parcellation (and network weight computation is skipped).
        Default: 'median'. Supported summary statistics include 'median', 'mode', 'mean', 'std', 'gmean', 'hmean', 'variation',
        'entropy', 'skew' and 'kurtosis'.

        Other appropriate summary statistics listed under scipy.stats could used
        by passing in a callable with their parameters encapsulated:
        https://docs.scipy.org/doc/scipy/reference/stats.html#statistical-functions
        For example, if you would like to compute 3rd k-statistic, you could construct a callable and passing ``third_kstat`` as in the argument:

        .. code-block:: python

            third_kstat  = lambda array: scipy.stats.kstat(array, n = 3)
            roi_medians = roiwise_stats_indiv(subject_id_list, fs_dir, base_feature, chosen_measure = third_kstat,
                atlas, fwhm, out_dir=None, return_results=True)

        Other possible options could trimmed mean estimator with 5% outliers removed or 3rd k-statistic:
        .. code-block:: python
            trimmed_mean = lambda array: scipy.stats.trim_mean(array, proportiontocut = 0.05)
            third_kstat  = lambda array: scipy.stats.kstat(array, n = 3)

        Notes: 'hmean' requires all values be positive.

    atlas : str
        Name of the atlas whose parcellation to be used.
        Available choices for cortical parcellation: ['FSAVERAGE', 'GLASSER2016'].
        Volumetric whole-brain atlases will be added soon.

    smoothing_param : scalar
        Smoothing parameter, which could be fwhm for Freesurfer cortical features,
        or another relevant for the chosen base_feature.
        Default: assumed as fwhm=10mm for the default feature choice 'thickness'

    node_size : scalar, optional
        Parameter to indicate the size of the ROIs, subparcels or patches, depending on type of atlas or feature.
        Not implemented.

    out_dir : str, optional
        Path to output directory to store results.
        Default: None, results are returned, but not saved to disk.
        If this is None, return_results must be true.

    return_results : bool
        Flag to indicating whether to keep the results to be returned to caller method.
        Helps to save memory (as it doesn't retain results all subjects and weight combinations),
        when running from command line interface (or HPC). Default: False
        If this is False, out_dir must be specified to save the results to disk.

    Returns
    -------
    roi_stats_all : dict, None
        If return_results is True, this will be a dictionary keyed in by subject_ID
        The value of each key roi_summary_all[subject] is
        a numpy array of length k, with k = number of nodes in the atlas parcellation.
        If return_results is False, this will be None, which is the default.
    """

    check_params_single_edge(base_feature, input_dir, atlas, smoothing_param, node_size, out_dir, return_results)
    subject_id_list, num_subjects, max_id_width, nd_id = check_subjects(subject_id_list)
    stat_func_list, stat_func_names, num_stats, max_stat_width, nd_st = check_stat_methods(chosen_roi_stats)

    # roi_labels, ctx_annot = parcellate.freesurfer_roi_labels(atlas)
    # uniq_rois, roi_size, num_nodes = roi_info(roi_labels)
    uniq_rois, centroids, roi_labels = parcellate.roi_labels_centroids(atlas)

    print('\nProcessing {} features resampled to {} atlas,'
          ' smoothed at {} with node size {}'.format(base_feature, atlas, smoothing_param, node_size))

    if return_results:
        roi_stats_all = dict()
    else:
        roi_stats_all = None
        if out_dir is None:
            raise ValueError('When return_results=False, out_dir must be specified to be able to save the results.')
        if not pexists(out_dir):
            os.mkdir(out_dir)

    for sub_idx, subject in enumerate(subject_id_list):

        try:
            features = import_features(input_dir, [subject, ], base_feature)
        except:
            raise IOError('Unable to read {} features for {}\n Skipping it.'.format(base_feature, subject))
            continue

        data, rois = mask_background_roi(features[subject], roi_labels, parcellate.null_roi_name)
        for ss, stat_func in enumerate(stat_func_list):
            sys.stdout.write('\nProcessing id {sid:{id_width}} ({sidnum:{nd_id}}/{numsub:{nd_id}}) -- '
                             'statistic {stname:{stat_name_width}} ({statnum:{nd_st}}/{numst:{nd_st}})'
                             ' :'.format(sid=subject, sidnum=sub_idx + 1, numsub=num_subjects,
                                         stname=stat_func_names[ss], statnum=ss + 1, numst=num_stats,
                                         id_width=max_id_width, stat_name_width=max_stat_width, nd_id=nd_id,
                                         nd_st=nd_st))

            try:
                roi_stats = calc_roi_statistics(data, rois, uniq_rois, stat_func)
                expt_id_no_network = stamp_experiment(base_feature, stat_func_names[ss], atlas, smoothing_param,
                                                      node_size)
                save_summary_stats(roi_stats, out_dir, subject, expt_id_no_network)
                sys.stdout.write('Done.')
            except KeyboardInterrupt:
                print('Exiting on keyborad interrupt! \n'
                      'Abandoning the remaining processing for {} stats:\n'
                      '{}.'.format(num_stats - ss, stat_func_names[ss:]))
                sys.exit(1)
            except:
                traceback.print_exc()
                logging.debug(
                    'Error : unable to compute roi-wise {} for {}. Skipping it.'.format(stat_func_names[ss], subject))

        if return_results:
            roi_stats_all[subject] = roi_stats

    return roi_stats_all


def calc_roi_statistics(data, rois, uniq_rois, given_callable=np.median):
    "Returns the requested ROI statistics."

    roi_stats = np.array([given_callable(data[rois == roi]) for roi in uniq_rois])

    return roi_stats


def import_features(input_dir,
                    subject_id_list,
                    base_feature,
                    fwhm=cfg.default_smoothing_param,
                    atlas=cfg.default_atlas):
    "Wrapper to support input data of multiple types and multiple packages."

    if isinstance(subject_id_list, str):
        subject_id_list = [subject_id_list, ]

    base_feature = base_feature.lower()
    if base_feature in cfg.features_freesurfer:
        features = freesurfer.import_features(input_dir, subject_id_list,
                                              base_feature=base_feature,
                                              fwhm=fwhm, atlas=atlas)
    elif base_feature in cfg.features_fsl:
        features = fsl_import(input_dir, subject_id_list, base_feature, fwhm=fwhm, atlas=atlas)
    else:
        raise NotImplementedError('Invalid or choice not implemented!\n'
                                  'Choose one of \n {}'.format(cfg.base_feature_list))

    return features


def fsl_import(input_dir,
               subject_id_list,
               base_feature,
               fwhm=cfg.default_smoothing_param,
               atlas=cfg.default_atlas):
    "To be implemented."

    if base_feature not in cfg.features_fsl:
        raise NotImplementedError

    return


def get_triu_handle_inf_nan(weights_matrix):
    "Issue a warning when NaNs or Inf are found."

    if weights_matrix is None:
        raise ValueError('Computation failed.')

    upper_tri_vec = weights_matrix[np.triu_indices_from(weights_matrix, 1)]

    warn_nan(upper_tri_vec)

    return upper_tri_vec


def mask_background_roi(data, labels, ignore_label):
    "Returns everything but specified label"

    if data.size != labels.size or data.shape != labels.shape:
        raise ValueError('features and membership/group labels differ in length/shape!')

    mask = labels != ignore_label
    masked_data    = data[mask]
    masked_labels  = labels[mask]

    if masked_data.size != masked_labels.size or masked_data.shape != masked_labels.shape:
        raise ValueError('features and membership/group labels (after removing background ROI) differ in length/shape!')

    return masked_data, masked_labels


def roi_info(roi_labels):
    "Unique ROIs in a given atlas parcellation, count and size. Excludes the background"

    uniq_rois_temp, roi_size_temp = np.unique(roi_labels, return_counts=True)

    # removing the background label
    index_bkgnd = np.argwhere(uniq_rois_temp == parcellate.null_roi_name)[0]
    uniq_rois = np.delete(uniq_rois_temp, index_bkgnd)
    roi_size = np.delete(roi_size_temp, index_bkgnd)

    num_nodes = len(uniq_rois)

    return uniq_rois, roi_size, num_nodes


def save_summary_stats(data_vec, out_dir, subject, str_suffix=None):
    "Saves the ROI medians to disk."

    if out_dir is not None:
        # get outpath returned from hiwenet, based on dist name and all other parameters
        # choose out_dir name  based on dist name and all other parameters
        out_subject_dir = pjoin(out_dir, subject)
        if not pexists(out_subject_dir):
            os.mkdir(out_subject_dir)

        if str_suffix is not None:
            out_file_name = '{}_roi_stats.csv'.format(str_suffix)
        else:
            out_file_name = 'roi_stats.csv'

        out_weights_path = pjoin(out_subject_dir, out_file_name)

        try:
            np.savetxt(out_weights_path, data_vec, fmt='%.5f')
            print('\nSaved roi stats to \n{}'.format(out_weights_path))
        except:
            print('\nUnable to save extracted features to {}'.format(out_weights_path))
            traceback.print_exc()

    return


def save_graph(graph_nx, out_dir, subject, str_suffix=None):
    "Saves the features to disk."

    if out_dir is not None:
        # get outpath returned from hiwenet, based on dist name and all other parameters
        # choose out_dir name  based on dist name and all other parameters
        out_subject_dir = pjoin(out_dir, subject)
        if not pexists(out_subject_dir):
            os.mkdir(out_subject_dir)

        if str_suffix is not None:
            out_file_name = '{}_graynet.graphml'.format(str_suffix)
        else:
            out_file_name = 'graynet.graphml'

        out_weights_path = pjoin(out_subject_dir, out_file_name)

        try:
            nx.info(graph_nx)
            nx.write_graphml(graph_nx, out_weights_path, encoding='utf-8')
            print('\nSaved the graph to \n{}'.format(out_weights_path))
        except:
            print('\nUnable to save graph to \n{}'.format(out_weights_path))
            traceback.print_exc()

    return


def save(weight_vec, out_dir, subject, str_suffix=None):
    "Saves the features to disk."

    if out_dir is not None:
        # get outpath returned from hiwenet, based on dist name and all other parameters
        # choose out_dir name  based on dist name and all other parameters
        out_subject_dir = pjoin(out_dir, subject)
        if not pexists(out_subject_dir):
            os.mkdir(out_subject_dir)

        if str_suffix is not None:
            out_file_name = '{}_graynet.csv'.format(str_suffix)
        else:
            out_file_name = 'graynet.csv'

        out_weights_path = pjoin(out_subject_dir, out_file_name)

        try:
            np.savetxt(out_weights_path, weight_vec, fmt='%.5f')
            print('\nSaved the features to \n{}'.format(out_weights_path))
        except:
            print('\nUnable to save features to {}'.format(out_weights_path))
            traceback.print_exc()

    return


def stamp_experiment(base_feature, method_name, atlas, smoothing_param, node_size):
    "Constructs a string to uniquely identify a given feature extraction method."

    # expt_id = 'feature_{}_atlas_{}_smoothing_{}_size_{}'.format(base_feature, atlas, smoothing_param, node_size)
    expt_id = '{}_{}_{}_smoothing{}_size{}'.format(method_name, base_feature, atlas, smoothing_param, node_size)

    return expt_id


def stamp_expt_weight(base_feature, atlas, smoothing_param, node_size, weight_method):
    "Constructs a string to uniquely identify a given feature extraction method."

    # expt_id = 'feature_{}_atlas_{}_smoothing_{}_size_{}'.format(base_feature, atlas, smoothing_param, node_size)
    expt_id = '{}_{}_smoothing{}_size{}_edgeweight_{}'.format(base_feature, atlas, smoothing_param, node_size,
                                                              weight_method)

    return expt_id


def cli_run():
    "command line interface!"

    subject_ids_path, input_dir, base_feature_list, \
    weight_method, do_multi_edge, summary_stat, num_bins, edge_range, \
    atlas, out_dir, node_size, smoothing_param, roi_stats, num_procs = parse_args()

    # when run from CLI, results will not be received
    # so no point in wasting memory maintaining a very big array
    return_results = False

    if do_multi_edge:
        from graynet.multi_edge import extract_multiedge
        print('Computing multiple edges ... ')
        extract_multiedge(subject_ids_path, input_dir,
                    base_feature_list=base_feature_list,
                    weight_method_list=weight_method,
                    summary_stat=summary_stat,
                    num_bins=num_bins, edge_range=edge_range,
                    atlas=atlas, smoothing_param=smoothing_param,
                    node_size=node_size, out_dir=out_dir,
                    return_results=return_results, num_procs=num_procs)
    else:
        base_feature = base_feature_list[0]
        if weight_method is not None:
            print('Computing single edge ... ')
            extract(subject_ids_path, input_dir,
                    base_feature=base_feature,
                    weight_method_list=weight_method,
                    num_bins=num_bins, edge_range=edge_range,
                    atlas=atlas, smoothing_param=smoothing_param,
                    node_size=node_size, out_dir=out_dir,
                    return_results=return_results, num_procs=num_procs)
        else:
            print('Computing ROI summary stats -- skipping computation of any network weights.')
            roiwise_stats_indiv(subject_ids_path, input_dir, base_feature,
                                roi_stats, atlas, smoothing_param, node_size,
                                out_dir, return_results)

    return


def get_parser():
    "Method to specify arguments and defaults. "

    help_text_subject_ids   = "Path to file containing list of subject IDs (one per line)"
    help_text_input_dir     = "Path to a folder containing input data. It could ,for example, be a Freesurfer SUBJECTS_DIR, if the chosen feature is from Freesurfer output."
    help_text_feature       = "Type of feature to be used for analysis. Default: '{}'. Choices: {}".format(cfg.default_feature_single_edge[0], cfg.base_feature_list)
    help_text_multi_edge    = "Option to compute multiple edges between ROIs based on different features. Default False. If True, two valid features must be specified."
    help_text_summary_stat  = "Summary statistic to compute of the weights from multiple edges." \
                              "This must be a string representing a method (like 'median', 'prod' or 'max'). " \
                              "that is available as a member of numpy or scipy.stats."
    help_text_weight        = "List of methods used to estimate the weight of the edge between the pair of nodes."  # .format(cfg.default_weight_method)
    help_text_num_bins      = "Number of bins used to construct the histogram within each ROI or group. Default : {}".format(cfg.default_num_bins)
    help_text_edge_range    = "The range of edges (two finite values) within which to bin the given values " \
                              "e.g. --edge_range 0.0 5.0 - this can be helpful (and important) to ensure correspondence " \
                              "across multiple invocations of graynet (for different subjects), " \
                              "in terms of range across all bins as well as individual bin edges. " \
                              "Default : {}, to automatically compute from the given values.".format(cfg.default_edge_range)

    help_text_roi_stats = "Option to compute summary statistics within each ROI of the chosen parcellation. " \
                          "These statistics (such as the median) can serve as a baseline for network-level values produced by graynet. " \
                          "Options for summary statistics include 'median', 'entropy', 'kurtosis' and " \
                          "any other appropriate summary statistics listed under scipy.stats: https://docs.scipy.org/doc/scipy/reference/stats.html#statistical-functions . " \
                          "When this option is chosen, network computation is not allowed. You need to compute networks and ROI stats separately."
    help_text_atlas     = "Name of the atlas to define parcellation of nodes/ROIs. Default: '{}'".format(cfg.default_atlas)
    help_text_parc_size = "Size of individual node for the atlas parcellation. Default : {}".format(cfg.default_node_size)
    help_text_smoothing = "Smoothing parameter for feature. Default: FWHM of {} for Freesurfer thickness".format(
        cfg.default_smoothing_param)

    help_text_num_procs = "Number of CPUs to use in parallel to speed up processing. Default : {}, capping at available number of CPUs in the processing node.".format(cfg.default_num_procs)

    parser = argparse.ArgumentParser(prog="graynet")

    parser.add_argument("-s", "--subject_ids_path", action="store", dest="subject_ids_path",
                        required=True,
                        help=help_text_subject_ids)

    parser.add_argument("-i", "--input_dir", action="store", dest="input_dir",
                        required=True, help=help_text_input_dir)

    parser.add_argument("-f", "--feature", action="store",
                        dest="features",
                        nargs='*',
                        default=cfg.default_feature_single_edge, required=False,
                        help=help_text_feature)

    parser.add_argument("-o", "--out_dir", action="store", dest="out_dir",
                        default=None, required=False,
                        help="Where to save the extracted features. ")

    method_selector = parser.add_argument_group(title='Type of computation',
                                                description='Choose one among single edge, multiedge or simply ROI stats.')
    # method_selector = parser.add_argument_group(required=True)
    method_selector.add_argument("-w", "--weight_method", action="store",
                                 dest="weight_methods",
                                 nargs='*',
                                 default=None, required=False, help=help_text_weight)

    method_selector.add_argument("-r", "--roi_stats", action="store", dest="roi_stats",
                                 nargs='*', default=None, help=help_text_roi_stats)

    method_selector.add_argument("-m", "--do_multi_edge", action="store_true", dest="do_multi_edge",
                                 default=False, required=False,
                                 help=help_text_multi_edge)

    method_selector.add_argument("-t", "--summary_stat", action="store", dest="summary_stat",
                                 default=cfg.multi_edge_summary_func_default, required=False,
                                 help=help_text_summary_stat)

    method_params = parser.add_argument_group(title='Weight parameters',
                                              description='Parameters relevant to histogram edge weight calculations')
    method_params.add_argument("-e", "--edge_range", action="store", dest="edge_range",
                               default=cfg.default_edge_range, required=False,
                               nargs=2, metavar=('min', 'max'),
                               help=help_text_edge_range)
    method_params.add_argument("-b", "--num_bins", action="store", dest="num_bins",
                                default=cfg.default_num_bins, required=False,
                                help=help_text_num_bins)


    atlas_params = parser.add_argument_group(title='Atlas',
                                             description="Parameters describing the atlas, its parcellation and any smoothing of features.")
    atlas_params.add_argument("-a", "--atlas", action="store", dest="atlas",
                              default=cfg.default_atlas, required=False,
                              help=help_text_atlas)

    atlas_params.add_argument("-n", "--node_size", action="store", dest="node_size",
                              default=cfg.default_node_size, required=False,
                              help=help_text_parc_size)

    atlas_params.add_argument("-p", "--smoothing_param", action="store", dest="smoothing_param",
                              default=cfg.default_smoothing_param, required=False,
                              help=help_text_smoothing)

    computing_params = parser.add_argument_group(title='Computing', description='Options related to computing and parallelization.')
    computing_params.add_argument('-c', '--num_procs', action='store', dest='num_procs',
                                  default=cfg.default_num_procs, required=False, help=help_text_num_procs)

    return parser


def parse_args():
    """Parser/validator for the cmd line args."""

    parser = get_parser()

    if len(sys.argv) < 2:
        parser.print_help()
        logging.warning('Too few arguments!')
        parser.exit(1)

    # parsing
    try:
        params = parser.parse_args()
    except Exception as exc:
        print(exc)
        raise ValueError('Unable to parse command-line arguments.')

    subject_ids_path = os.path.abspath(params.subject_ids_path)
    if not os.path.exists(subject_ids_path):
        raise IOError("Given subject IDs file doesn't exist.")

    input_dir = os.path.abspath(params.input_dir)
    if not os.path.exists(input_dir):
        raise IOError("Given input directory doesn't exist.")

    out_dir = params.out_dir
    if out_dir is not None:
        if not pexists(out_dir):
            os.mkdir(out_dir)

    feature_list = utils.check_features(params.features)

    do_multi_edge = bool(params.do_multi_edge)
    summary_stat = params.summary_stat
    if do_multi_edge:
        # ensure atleast two features
        if len(feature_list) < 2:
            raise ValueError('To enable multi-edge computation, specify atleast two valid features.')

        utils.check_stat_methods(summary_stat)
    else:
        summary_stat = None
        if len(feature_list) > 1:
            raise ValueError('For single edge computation, only one feature can be specified.')

    # validating choices and doing only one of the two
    weight_methods = params.weight_methods
    roi_stats = params.roi_stats
    if weight_methods is not None:
        weight_method_list, _, _, _ = check_weights(weight_methods)
        if roi_stats is not None:
            print('ROI stats requested with network weights computation - not allowed.')
            sys.exit(1)
        roi_stats = None
    elif roi_stats is not None:
        roi_stats, _, _, _, _ = check_stat_methods(roi_stats)
        weight_method_list = None
    else:
        raise ValueError('One of weight_method and roi_stats must be chosen.')

    # num_procs will be validated inside in the functions using it.

    return subject_ids_path, input_dir, \
           feature_list, weight_method_list, \
           do_multi_edge, summary_stat, \
           params.num_bins, params.edge_range, \
           params.atlas, out_dir, params.node_size, params.smoothing_param, roi_stats, params.num_procs


if __name__ == '__main__':
    cli_run()
