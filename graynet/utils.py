"""

Common utilities

"""
import collections
import logging
import os
import sys
import traceback
from genericpath import exists as pexists
from os.path import join as pjoin, isfile, realpath, getsize
import nibabel
import networkx as nx
from multiprocessing import cpu_count

import numpy as np

from graynet import config_graynet as cfg, freesurfer


def check_features(base_feature_list):
    "Validates the choice of features"

    if base_feature_list in [None, '']:
        raise ValueError('feature list can not be empty.')

    # when a string is specified, making it a list
    if isinstance(base_feature_list, str):
        base_feature_list = [base_feature_list, ]

    given_list = unique_order(base_feature_list)
    given = set(given_list)
    allowed = set(cfg.base_feature_list)
    if not given.issubset(allowed):
        unrecog_methods = given.difference(allowed)
        raise NotImplementedError('features unrecognized: \n {}\n'
                                  ' choose one of :\n\t {}'
                                  ''.format(unrecog_methods, allowed))

    return given_list


def check_atlas(atlas):
    """Validation of atlas input."""

    # when its a name for pre-defined atlas
    if isinstance(atlas, str):
        if not pexists(atlas):  # just a name
            atlas = atlas.lower()
            if atlas not in cfg.atlas_list:
                raise ValueError(
                    'Invalid choice of atlas {}.'
                    ' Accepted : {}'.format(atlas, cfg.atlas_list))
        elif os.path.isdir(atlas):  # cortical atlas in Freesurfer org
            if not check_atlas_annot_exist(atlas):
                raise ValueError(
                    'Given atlas folder does not contain Freesurfer label annot files. '
                    'Needed : given_atlas_dir/label/?h.aparc.annot')
        elif pexists(atlas):  # may be a volumetric atlas?
            try:
                atlas = nibabel.load(atlas)
            except:
                traceback.print_exc()
                raise ValueError('Unable to read the provided image volume. '
                                 'Must be a nifti 2d volume, readable by nibabel.')
        else:
            raise ValueError('Unable to decipher or use the given atlas.')
    else:
        raise NotImplementedError('Atlas must be a string, providing a name or '
                                  'path to Freesurfer folder or a 3D nifti volume.')

    return atlas


def unique_order(seq):
    "Removes duplicates while preserving order"

    uniq = list()
    for element in seq:
        if element not in uniq:
            uniq.append(element)

    return uniq


def make_output_path_graph(out_dir, subject, str_prefixes):
    "Constructs path to save a multigraph to disk."

    if out_dir is not None:
        # get outpath returned from hiwenet, based on dist name and all other params
        # choose out_dir name  based on dist name and all other parameters
        out_subject_dir = pjoin(out_dir, subject)
        if not pexists(out_subject_dir):
            os.mkdir(out_subject_dir)

        if isinstance(str_prefixes, str):
            str_prefixes = [str_prefixes, ]

        out_file_name = '{}_graynet.graphml'.format('_'.join(str_prefixes))
        out_weights_path = pjoin(out_subject_dir, out_file_name)
    else:
        out_weights_path = None

    return out_weights_path


def save_graph(graph, out_path, identifier=''):
    "Saves the given graph to disk."

    if out_path is not None:
        try:
            nx.write_graphml(graph, out_path, encoding='utf-8')
            print('\nSaved the {} graph to \n{}'.format(identifier, out_path))
        except:
            print('\nUnable to save {} graph to \n{}'.format(identifier, out_path))
            traceback.print_exc()

    return out_path


def check_num_procs(num_procs=cfg.default_num_procs):
    "Ensures num_procs is finite and <= available cpu count."

    num_procs = int(num_procs)
    avail_cpu_count = cpu_count()
    if num_procs < 1 or not np.isfinite(num_procs) or num_procs is None:
        num_procs = 1
        print('Invalid value for num_procs. Using num_procs=1')
    elif num_procs > avail_cpu_count:
        print('# CPUs requested higher than available '
              '- choosing {}'.format(avail_cpu_count))
        num_procs = avail_cpu_count

    return num_procs


def check_stat_methods(stat_list=None):
    "Validates the choice and returns a callable to compute summaries."

    from scipy import stats as sp_stats
    from functools import partial

    if stat_list is None:
        stat_list = [np.median, ]

    if not isinstance(stat_list, (list, tuple, set)):
        # when a single method is specified by a str or callable
        if isinstance(stat_list, str) or callable(stat_list):
            stat_list = [stat_list, ]
        else:
            raise ValueError('Unrecognized stat method: must be a str or callable '
                             'or a list of them.')

    stat_callable_list = list()
    for stat in stat_list:
        if isinstance(stat, str):
            stat = stat.lower()
            if hasattr(np, stat):
                summary_callable = getattr(np, stat)
            elif hasattr(sp_stats, stat):
                summary_callable = getattr(sp_stats, stat)
            else:
                raise AttributeError('Chosen measure {} is not a member of numpy '
                                     'or scipy.stats.'.format(stat))
        elif callable(stat):
            summary_callable = stat
        else:
            raise ValueError('summary measure is not recognized.')

        stat_callable_list.append(summary_callable)

    # constructing names
    names_callable = list()
    for func in stat_callable_list:
        try:
            method_name = func.__name__
        except:
            if isinstance(func, partial):
                method_name = func.func.__name__
            else:
                raise ValueError('name of callable {} could not be obtained'.format(func))
        method_name = method_name.replace(' ', '_')
        names_callable.append(method_name)

    num_stats = len(stat_callable_list)
    num_digits_stat_size = len(str(num_stats))
    max_wtname_width = max(map(len, names_callable))

    return stat_callable_list, names_callable, num_stats, max_wtname_width, num_digits_stat_size


def warn_nan(array):
    "Raises a warning when non-finite or NaN values are found."

    num_nonfinite = np.count_nonzero(np.logical_not(np.isfinite(array)))
    if num_nonfinite > 0:
        logging.warning('{} non-finite values are found.'.format(num_nonfinite))

    return


def check_subjects(subjects_info):
    "Ensure subjects are provided and their data exist."

    if isinstance(subjects_info, str):
        if not pexists(subjects_info):
            raise IOError('path to subject list does not exist: {}'.format(subjects_info))
        subjects_list = np.genfromtxt(subjects_info, dtype=str)
    elif isinstance(subjects_info, collections.Iterable):
        if len(subjects_info) < 1:
            raise ValueError('Empty subject list.')
        subjects_list = subjects_info
    else:
        raise ValueError('Invalid value provided for subject list. \n '
                         'Must be a list of paths, or '
                         'path to a file containing one path per line for each subject.')

    subject_id_list = np.atleast_1d(subjects_list)
    num_subjects = subject_id_list.size
    if num_subjects < 1:
        raise ValueError('Input subject list is empty.')

    num_digits_id_size = len(str(num_subjects))
    max_id_width = max(map(len, subject_id_list))

    return subject_id_list, num_subjects, max_id_width, num_digits_id_size


def check_weights(weight_method_list):
    "Ensures weights are implemented and atleast one choice is given."

    if isinstance(weight_method_list, str):
        weight_method_list = [weight_method_list, ]

    if isinstance(weight_method_list, collections.Iterable):
        if len(weight_method_list) < 1:
            raise ValueError('Empty weight list. Atleast one weight must be provided.')
    else:
        raise ValueError('Weights list must be an iterable. Given: {}'
                         ''.format(type(weight_method_list)))

    for weight in weight_method_list:
        if weight not in cfg.implemented_weights:
            raise NotImplementedError('Method {} not implemented. '
                                      'Choose one of : \n {}'
                                      ''.format(weight, cfg.implemented_weights))

    num_weights = len(weight_method_list)
    num_digits_wm_size = len(str(num_weights))
    max_wtname_width = max(map(len, weight_method_list))

    return weight_method_list, num_weights, max_wtname_width, num_digits_wm_size


def check_edge_range(edge_range_spec):
    "Validates the edge rage specified"

    if edge_range_spec is None:
        edge_range = edge_range_spec
    elif isinstance(edge_range_spec, (collections.Sequence, np.ndarray)):
        if len(edge_range_spec) != 2:
            raise ValueError('edge_range must be a tuple of two values: (min, max)')
        if edge_range_spec[0] >= edge_range_spec[1]:
            raise ValueError('edge_range : min {} is not less than max {} !'
                             ''.format(edge_range_spec[0], edge_range_spec[1]))

        # CLI args are strings unless converted to numeric
        edge_range = np.float64(edge_range_spec)
        if not np.all(np.isfinite(edge_range)):
            raise ValueError('Infinite or NaN values in edge range : {}'.format(edge_range_spec))

        # converting it to tuple to make it immutable
        edge_range = tuple(edge_range)
    else:
        raise ValueError('Invalid edge range! '
                         'Must be a tuple of two values (min, max)')

    return edge_range


def check_num_bins(num_bins):
    "Validates the number of bins chosen"

    if isinstance(num_bins, str):
        # possible when called from CLI
        num_bins = np.float(num_bins)

    # rounding it to ensure it is int
    num_bins = np.rint(num_bins)

    if np.isnan(num_bins) or np.isinf(num_bins):
        raise ValueError('Invalid value for number of bins! '
                         'Choose a natural number >= {}'
                         ''.format(cfg.default_minimum_num_bins))

    return num_bins


def check_weight_params(num_bins, edge_range_spec):
    "Ensures parameters are valid and type casts them."

    num_bins = check_num_bins(num_bins)

    edge_range = check_edge_range(edge_range_spec)

    return num_bins, edge_range


def check_edge_range_dict(edge_range_dict, base_feature_list,
                          predefined_ranges=cfg.edge_range_predefined):
    "Ensures ranges were specified for each feature, and they are valid or automatic(None)"

    print('Setting given edge range ...')
    for feature in base_feature_list:
        sys.stdout.write('\n <---- {:20s} '.format(feature))
        if feature in edge_range_dict:
            edge_range_dict[feature] = check_edge_range(edge_range_dict[feature])
            sys.stdout.write(': {} ----> \n'.format(edge_range_dict[feature]))
        elif feature in predefined_ranges:
            sys.stdout.write(
                'edge range not given! Using predefined: {}'.format(predefined_ranges[feature]))
            edge_range_dict[feature] = predefined_ranges[feature]
        else:
            # covers the case of edge_range_dict being None
            sys.stdout.write('edge range not given or predefined! '
                             'Setting it automatic (may change for each subject)')
            edge_range_dict[feature] = None

    return edge_range_dict


def check_params_single_edge(base_features, in_dir, atlas, smoothing_param,
                             node_size, out_dir, return_results):
    """"""

    check_features(base_features)

    check_atlas(atlas)

    if not pexists(in_dir):
        raise IOError('Input directory at {} does not exist.'.format(in_dir))

    if out_dir is None and return_results is False:
        raise ValueError('Results are neither saved to disk, '
                         'nor being received when returned!\n'
                         'Specify out_dir (not None) or make return_results=True')

    if out_dir is not None and not pexists(out_dir):
        os.mkdir(out_dir)

    # no checks on subdivison size yet, as its not implemented

    return


def stamp_expt_multiedge(base_feature_list, atlas, smoothing_param, node_size,
                         weight_method):
    "Constructs a string to uniquely identify a given experiment."

    import re
    all_words = re.split('_|; |, |\*|\n| ', ' '.join(base_feature_list))
    feat_repr = '_'.join(unique_order(all_words))
    expt_id = '{}_{}_smth{}_{}_{}'.format(feat_repr, atlas, smoothing_param, node_size,
                                          weight_method)

    return expt_id


def check_params_multiedge(base_feature_list, input_dir, atlas, smoothing_param,
                           node_size, out_dir, return_results):
    """Validation of parameters and appropriate type casting if necessary."""

    check_features(base_feature_list)

    check_atlas(atlas)

    if not pexists(input_dir):
        raise IOError('Input directory at {} does not exist.'.format(input_dir))

    if out_dir is None and return_results is False:
        raise ValueError('Results are neither saved to disk, '
                         'nor being received when returned.\n'
                         'Specify out_dir (not None) or make return_results=True')

    if out_dir is not None and not pexists(out_dir):
        os.mkdir(out_dir)

    # no checks on subdivison size yet, as its not implemented

    return


def calc_roi_statistics(data, rois, uniq_rois, given_callable=np.median):
    "Returns the requested ROI statistics."

    roi_stats = np.array([given_callable(data[rois == roi]) for roi in uniq_rois])

    return roi_stats


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
        raise ValueError('features and membership (group labels) differ'
                         ' in length or shape!')

    mask = labels != ignore_label
    masked_data = data[mask]
    masked_labels = labels[mask]

    if masked_data.size != masked_labels.size or \
            masked_data.shape != masked_labels.shape:
        raise ValueError('features and membership (group labels), '
                         'after removing background ROI, differ in length/shape!')

    return masked_data, masked_labels


def roi_info(roi_labels, freesurfer_annot=True):
    "Unique ROIs in a given atlas parcellation, count and size. Excludes the background"

    uniq_rois_temp, roi_size_temp = np.unique(roi_labels, return_counts=True)

    # removing the background label
    if freesurfer_annot:
        index_bkgnd = np.argwhere(uniq_rois_temp == cfg.null_roi_name)[0]
    else:
        index_bkgnd = np.argwhere(uniq_rois_temp == cfg.null_roi_index)[0]

    uniq_rois = np.delete(uniq_rois_temp, index_bkgnd)
    roi_size = np.delete(roi_size_temp, index_bkgnd)

    num_nodes = len(uniq_rois)

    return uniq_rois, roi_size, num_nodes


def check_atlas_annot_exist(atlas_dir, hemi_list=None):
    " Checks for the presence of atlas annotations "

    if hemi_list is None:
        hemi_list = ['lh', 'rh']

    for hemi in hemi_list:
        annot_path = pjoin(atlas_dir, 'label', '{}.aparc.annot'.format(hemi))
        if not pexists(annot_path) or os.path.getsize(annot_path) == 0:
            return False

    return True


def stamp_experiment(base_feature, method_name, atlas, smoothing_param, node_size):
    "Constructs a string to uniquely identify a given feature extraction method."

    # expt_id = 'feature_{}_atlas_{}_smoothing_{}_size_{}'.format(base_feature, atlas, smoothing_param, node_size)
    expt_id = '{}_{}_{}_smoothing{}_size{}'.format(method_name, base_feature, atlas,
                                                   smoothing_param, node_size)

    return expt_id


def stamp_expt_weight(base_feature, atlas, smoothing_param, node_size, weight_method):
    "Constructs a string to uniquely identify a given feature extraction method."

    # expt_id = 'feature_{}_atlas_{}_smoothing_{}_size_{}'.format(base_feature, atlas, smoothing_param, node_size)
    expt_id = '{}_{}_smoothing{}_size{}_edgeweight_{}' \
              ''.format(base_feature, atlas, smoothing_param, node_size,
                        weight_method)

    return expt_id


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
        features = fsl_import(input_dir, subject_id_list, base_feature,
                              fwhm=fwhm, atlas=atlas)
    elif base_feature in cfg.features_spm_cat:
        features = spm_cat_import(input_dir, subject_id_list, base_feature, atlas=atlas)
    else:
        raise NotImplementedError('Invalid or choice not implemented!\n'
                                  'Choose one of \n {}'.format(cfg.base_feature_list))

    return features


def save_summary_stats(roi_values, roi_labels, stat_name, out_dir, subject, str_suffix=None):
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
            with open(out_weights_path, 'w') as of:
                of.write('#roi,{}\n'.format(stat_name))
                for name, value in zip(roi_labels, roi_values):
                    of.write('{},{}\n'.format(name, value))
            # np.savetxt(out_weights_path, roi_values, fmt='%.5f')
            print('\nSaved roi stats to \n{}'.format(out_weights_path))
        except:
            print('\nUnable to save extracted features to {}'.format(out_weights_path))
            traceback.print_exc()

    return


def save_per_subject_graph(graph_nx, out_dir, subject, str_suffix=None):
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


def spm_cat_import(input_dir,
                   subject_id_list,
                   base_feature,
                   atlas=cfg.default_vbm_atlas):
    """Imports the voxelwise features for a given list of subjecct IDs."""

    features= dict()
    for sid in subject_id_list:
        try:
            print('Reading {} for {} ... '.format(base_feature, sid), end='')
            features[sid] = get_CAT_data(input_dir, sid, base_feature)
            print(' Done.')
        except:
            traceback.print_exc()
            raise ValueError('{} data for {} could not be read!'
                             ''.format(base_feature, sid))

    return features


def get_CAT_data(input_dir, sid, base_feature):
    """Returns the values in a specified image!"""

    img_path = get_SPM_CAT_img_path(input_dir, sid, base_feature)
    img = nibabel.load(img_path).get_data()

    return img

def get_SPM_CAT_img_path(input_dir, sid, base_feature):
    """Constructs the path for a given subject ID and feature"""

    return pjoin(input_dir, 'mri', '{}{}.nii'.format(
            cfg.features_spm_cat_prefixes[base_feature], sid))


def fsl_import(input_dir,
               subject_id_list,
               base_feature,
               fwhm=cfg.default_smoothing_param,
               atlas=cfg.default_atlas):
    "To be implemented."

    if base_feature not in cfg.features_fsl:
        raise NotImplementedError

    return