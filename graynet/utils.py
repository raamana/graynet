"""

Common utilities

"""
import collections
import logging
import os
import traceback
from genericpath import exists as pexists

import nibabel
from multiprocessing import cpu_count

import numpy as np

from graynet import config_graynet as cfg, parcellate


def check_features(base_feature_list):
    "Validates the choice of features"

    if base_feature_list in [None, '']:
        raise ValueError('feature list can not be empty.')

    # when a string is specified, making it a list
    if isinstance(base_feature_list,str):
        base_feature_list = [base_feature_list, ]

    # remove duplicates, preserve order
    given_list = []
    for feat in base_feature_list:
        if feat not in given_list:
            given_list.append(feat)

    given = set(given_list)
    allowed = set(cfg.base_feature_list)
    if not given.issubset(allowed):
        unrecog_methods = given.difference(allowed)
        raise NotImplementedError('Methods unrecognized: \n {}'.format(unrecog_methods))

    return given_list


def check_atlas(atlas):
    """Validation of atlas input."""

    # when its a name for pre-defined atlas
    if isinstance(atlas, str):
        if not pexists(atlas): # just a name
            atlas = atlas.lower()
            if atlas not in parcellate.atlas_list:
                raise ValueError('Invalid choice of atlas. Accepted : {}'.format(parcellate.atlas_list))
        elif os.path.isdir(atlas): # cortical atlas in Freesurfer org
            if not parcellate.check_atlas_annot_exist(atlas):
                raise ValueError('Given atlas folder does not contain Freesurfer label annot files. '
                                 'Needed : given_atlas_dir/label/?h.aparc.annot')
        elif pexists(atlas): # may be a volumetric atlas?
            try:
                atlas = nibabel.load(atlas)
            except:
                traceback.print_exc()
                raise ValueError('Unable to read the provided image volume. Must be a nifti 2d volume, readable by nibabel.')
        else:
            raise ValueError('Unable to decipher or use the given atlas.')
    else:
        raise NotImplementedError('Atlas must be a string, providing a name or path to Freesurfer folder or a 3D nifti volume.')

    return atlas


def check_num_procs(num_procs=cfg.default_num_procs):
    "Ensures num_procs is finite and <= available cpu count."

    num_procs  = int(num_procs)
    avail_cpu_count = cpu_count()
    if num_procs < 1 or not np.isfinite(num_procs) or num_procs is None:
        num_procs = 1
        print('Invalid value for num_procs. Using num_procs=1')
    elif num_procs > avail_cpu_count:
        print('# CPUs requested higher than available - choosing {}'.format(avail_cpu_count))
        num_procs = avail_cpu_count

    return num_procs


def check_stat_methods(stat_list=None):
    "Validates the choice and returns a callable to compute summaries."

    from scipy import stats as sp_stats
    from functools import partial

    if stat_list is None:
        stat_list = [np.median, ]

    if not isinstance(stat_list, list):
        # when a single method is specified by a str or callable
        if isinstance(stat_list, str) or callable(stat_list):
            stat_list = [stat_list, ]
        else:
            raise ValueError('Unrecognized stat method: must be a str or callable or a list of them.')

    stat_callable_list = list()
    for stat in stat_list:
        if isinstance(stat, str):
            stat = stat.lower()
            if hasattr(np, stat):
                summary_callable = getattr(np, stat)
            elif hasattr(sp_stats, stat):
                summary_callable = getattr(sp_stats, stat)
            else:
                raise AttributeError('Chosen measure {} is not a member of numpy or scipy.stats.'.format(stat))
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
                         'Must be a list of paths, or path to file containing list of paths, one for each subject.')

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
        raise ValueError('Weights list must be an iterable. Given: {}'.format(type(weight_method_list)))

    for weight in weight_method_list:
        if weight not in cfg.implemented_weights:
            raise NotImplementedError('Method {} not implemented. '
                                      'Choose one of : \n {}'.format(weight, cfg.implemented_weights))

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
            raise ValueError('edge_range : min {} is not less than max {} !'.format(edge_range_spec[0], edge_range_spec[1]))

        # CLI args are strings unless converted to numeric
        edge_range = np.float64(edge_range_spec)
        if not np.all(np.isfinite(edge_range)):
            raise ValueError('Infinite or NaN values in edge range : {}'.format(edge_range_spec))

        # converting it to tuple to make it immutable
        edge_range = tuple(edge_range)
    else:
        raise ValueError('Invalid edge range! Must be a tuple of two values (min, max)')

    return edge_range


def check_num_bins(num_bins):
    "Validates the number of bins chosen"

    if isinstance(num_bins, str):
        # possible when called from CLI
        num_bins = np.float(num_bins)

    # rounding it to ensure it is int
    num_bins = np.rint(num_bins)

    if np.isnan(num_bins) or np.isinf(num_bins):
        raise ValueError('Invalid value for number of bins! Choose a natural number >= {}'.format(cfg.default_minimum_num_bins))

    return num_bins


def check_weight_params(num_bins, edge_range_spec):
    "Ensures parameters are valid and type casts them."

    num_bins = check_num_bins(num_bins)

    edge_range = check_edge_range(edge_range_spec)

    return num_bins, edge_range


def check_edge_range_dict(edge_range_dict, base_feature_list, predefined_ranges=cfg.edge_range_predefined):
    "Ensures ranges were specified for each feature, and they are valid or automatic(None)"

    for feature in base_feature_list:
        print(' <---- {:20s} ----> '.format(feature))
        if feature in edge_range_dict:
            edge_range_dict[feature] = check_edge_range(edge_range_dict[feature])
            print('Setting given edge range : {}'.format(edge_range_dict[feature]))
        elif feature in predefined_ranges:
            print('edge range not given! Using predefined: {}'.format(predefined_ranges[feature]))
            edge_range_dict[feature] = predefined_ranges[feature]
        else:
            # covers the case of edge_range_dict being None
            print('edge range not given or predefined! Setting it automatic (may change for each subject)')
            edge_range_dict[feature] = None

    return edge_range_dict


def check_params_single_edge(base_features, in_dir, atlas, smoothing_param, node_size, out_dir, return_results):
    """"""

    check_features(base_features)

    if atlas.lower() not in parcellate.atlas_list:
        raise ValueError('Invalid atlas choice. Use one of {}'.format(parcellate.atlas_list))

    if not pexists(in_dir):
        raise IOError('Input directory at {} does not exist.'.format(in_dir))

    if out_dir is None and return_results is False:
        raise ValueError('Results are neither saved to disk or being received when returned.\n'
                         'Specify out_dir (not None) or make return_results=True')

    if out_dir is not None and not pexists(out_dir):
        os.mkdir(out_dir)

    # no checks on subdivison size yet, as its not implemented

    return


def stamp_expt_multiedge(base_feature_list, atlas, smoothing_param, node_size, weight_method):
    "Constructs a string to uniquely identify a given experiment."

    import re
    all_words = re.split('_|; |, |\*|\n| ', ' '.join(base_feature_list))
    feat_repr = '_'.join(set(all_words))
    expt_id   = '{}_{}_smth{}_{}_{}'.format(feat_repr, atlas, smoothing_param, node_size, weight_method)

    return expt_id


def check_params_multiedge(base_feature_list, input_dir, atlas, smoothing_param,
                           node_size, out_dir, return_results):
    """Validation of parameters and appropriate type casting if necessary."""

    check_features(base_feature_list)

    if atlas.lower() not in parcellate.atlas_list:
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