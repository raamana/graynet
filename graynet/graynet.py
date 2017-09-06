
__all__ = ['extract', 'roi_info', 'cli_run']

import collections
import os
import sys
import argparse
import traceback
import logging
from os.path import join as pjoin, exists as pexists

import hiwenet
import nibabel
import numpy as np

from sys import version_info

if version_info.major==2 and version_info.minor==7 and version_info.micro==13:
    import freesurfer
    import parcellate
elif version_info.major > 2:
    from . import freesurfer
    from . import parcellate
else:
    raise NotImplementedError('hiwenet supports only 2.7.13 or 3+. Upgrade to Python 3+ is recommended.')


np.seterr(divide='ignore', invalid='ignore')

__features_freesurfer = ['freesurfer_thickness', 'freesurfer_curv', 'freesurfer_sulc']
__features_fsl = ['gmdensity', ]

__base_feature_list = __features_freesurfer + __features_fsl

__default_weight_method = ( 'minowski', 'manhattan' )
__accepted_weight_list = [
    'chebyshev', 'chebyshev_neg', 'chi_square',
    'correlate', 'correlate_1',
    'cosine', 'cosine_1', 'cosine_2', 'cosine_alt',
    'euclidean', 'fidelity_based',
    'histogram_intersection', 'histogram_intersection_1',
    'jensen_shannon', 'kullback_leibler', 'manhattan', 'minowski',
    'noelle_1', 'noelle_2', 'noelle_3', 'noelle_4', 'noelle_5',
    'relative_bin_deviation', 'relative_deviation']
__default_num_bins = 25
__default_trim_percentile = 5

__default_feature = 'freesurfer_thickness'
__default_atlas = 'GLASSER2016'
__default_smoothing_param = 10
__default_node_size = None

def extract(subject_id_list, input_dir,
            base_feature = __default_feature,
            weight_method_list = __default_weight_method,
            atlas = __default_atlas, smoothing_param = __default_smoothing_param,
            node_size = __default_node_size, out_dir=None, return_results = False):
    """
    Extracts weighted networks from gray matters features based on Freesurfer processing.

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
    weight_method_list : list of str
        Names of covariance metrics to use to compute weights.
        Currently only those supported by hiwenet, which can be one of:
        [ 'chebyshev', 'chebyshev_neg', 'chi_square', 'correlate', 'correlate_1',
        'cosine', 'cosine_1', 'cosine_2', 'cosine_alt', 'euclidean', 'fidelity_based',
        'histogram_intersection', 'histogram_intersection_1', 'jensen_shannon', 'kullback_leibler',
        'manhattan', 'minowski', 'noelle_1', 'noelle_2', 'noelle_3', 'noelle_4', 'noelle_5',
        'relative_bin_deviation', 'relative_deviation']
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
        Not implemented.
    out_dir : str, optional
        Path to output directory to store results.
        Default: None, results are returned, but not saved to disk.
        If this is None, return_results must be true.
    return_results : bool
        Flag to indicating whether to keep the results to be returned to caller method.
        Helps to save memory (as it doesn't retain results all subjects and weight combinations),
        when running from commmand line interface (or HPC). Default: False
        If this is False, out_dir must be specified to save the results to disk.

    Returns
    -------
    edge_weights_all : dict, None
        If return_results is True, this will be a dictionary keyed in by a tuple: (weight method, subject_ID)
        The value of each edge_weights_all[(weight method, subject_ID)] is
        a numpy array of length p = k*(k-1)/2, with k = number of nodes in the atlas parcellation.
        If return_results is False, this will be None, which is the default.
    """

    __parameter_check(base_feature, input_dir, atlas, smoothing_param, node_size, out_dir, return_results)
    subject_id_list, num_subjects, max_id_width, nd_id = __subject_check(subject_id_list)
    weight_method_list, num_weights, max_wtname_width, nd_wm = __weight_check(weight_method_list)

    roi_labels, ctx_annot = parcellate.freesurfer_roi_labels(atlas)
    uniq_rois, roi_size, num_nodes = roi_info(roi_labels)

    print('\nProcessing {} features resampled to {} atlas,'
          ' smoothed at {} with node size {}'.format(base_feature, atlas, smoothing_param, node_size))

    if return_results:
        edge_weights_all = dict()
    else:
        edge_weights_all = None

    # edge_weights_all[weight_method] = dict() # np.zeros([num_subjects, np.int64(num_nodes*(num_nodes-1)/2)])
    for ss, subject in enumerate(subject_id_list):

        try:
            features = import_features(input_dir, [subject, ], base_feature)
        except:
            raise IOError('Unable to read {} features for {}\n Skipping it.'.format(base_feature, subject))
            continue

        data, rois = __remove_background_roi(features[subject], roi_labels, parcellate.null_roi_name)

        for ww, weight_method in enumerate(weight_method_list):
            # unique stamp for each subject and weight
            expt_id = __stamp_experiment(base_feature, atlas, smoothing_param, node_size, weight_method)
            sys.stdout.write('\nProcessing id {:{id_width}} ({:{nd_id}}/{:{nd_id}}) -- '
                             'weight {:{wtname_width}} ({:{nd_wm}}/{:{nd_wm}})'
                             ' :'.format(subject, ss+1, num_subjects, weight_method, ww+1, num_weights,
                                         id_width=max_id_width, wtname_width=max_wtname_width,
                                         nd_id=nd_id, nd_wm=nd_wm)) # controlling for number of digits

            # actual computation of pair-wise features
            try:
                edge_weights = hiwenet.extract(data, rois, weight_method)

            except (RuntimeError, RuntimeWarning) as runexc:
                print(runexc)
            except KeyboardInterrupt:
                print('Exiting on keyborad interrupt! \n'
                      'Abandoning the remaining processing for {} weights:\n'
                      '{}.'.format(num_weights-ww, weight_method_list[ww:]))
                sys.exit(1)
            except:
                print('Unable to extract {} features for {}'.format(weight_method, subject))
                traceback.print_exc()

            # saving to disk
            try:
                weight_vec = __get_triu_handle_inf_nan(edge_weights)
                __save(weight_vec, out_dir, subject, expt_id)
            except:
                raise IOError('Unable to save the vectorized features to:\n{}'.format(out_dir))

            # saving the results to memory only if needed.
            if return_results:
                edge_weights_all[(weight_method, subject)] = weight_vec

            sys.stdout.write('Done.')

    return edge_weights_all


def import_features(input_dir, subject_id_list, base_feature):
    "Wrapper to support input data of multiple types and multiple packages."

    if isinstance(subject_id_list, str):
        subject_id_list = [subject_id_list, ]

    base_feature = base_feature.lower()
    if base_feature in __features_freesurfer:
        features = freesurfer.import_features(input_dir, subject_id_list, base_feature)
    elif base_feature in __features_fsl:
        features = fsl_import(input_dir, subject_id_list, base_feature)
    else:
        raise NotImplementedError('Invalid or choice not implemented!\n'
                                  'Choose one of \n {}'.format(__base_feature_list))

    return features


def fsl_import(input_dir, subject_id_list, base_feature):
    "To be implemented."

    if base_feature not in __features_fsl:
        raise NotImplementedError

    return


def __get_triu_handle_inf_nan(weights_matrix):
    "Issue a warning when NaNs or Inf are found."

    if weights_matrix is None:
        raise ValueError('Computation failed.')

    upper_tri_vec = weights_matrix[np.triu_indices_from(weights_matrix, 1)]

    num_nonfinite = np.count_nonzero(np.logical_not(np.isfinite(upper_tri_vec)))
    if num_nonfinite > 0:
        logging.warning(' {} non-finite values are found.'.format(num_nonfinite))

    return upper_tri_vec


def __subject_check(subjects_info):
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
    num_digits_id_size = len(str(num_subjects))
    max_id_width = max(map(len, subject_id_list))

    return subject_id_list, num_subjects, max_id_width, num_digits_id_size


def __weight_check(weight_method_list):
    "Ensures weights are implemented and atleast one choice is given."

    if isinstance(weight_method_list, str):
        weight_method_list = [weight_method_list, ]

    if isinstance(weight_method_list, collections.Iterable):
        if len(weight_method_list) < 1:
            raise ValueError('Empty weight list. Atleast one weight must be provided.')
    else:
        raise ValueError('Weights list must be an iterable. Given: {}'.format(type(weight_method_list)))

    for weight in weight_method_list:
        if weight not in __accepted_weight_list:
            raise NotImplementedError('Method {} not implemented. '
                                      'Choose one of : \n {}'.format(weight, __accepted_weight_list))

    num_weights = len(weight_method_list)
    num_digits_wm_size = len(str(num_weights))
    max_wtname_width = max(map(len, weight_method_list))

    return weight_method_list, num_weights, max_wtname_width, num_digits_wm_size


def __remove_background_roi(data,labels, ignore_label):
    "Returns everything but specified label"

    mask = labels != ignore_label

    return data[mask], labels[mask]


def roi_info(roi_labels):
    "Unique ROIs in a given atlas parcellation, count and size. Excludes the background"

    uniq_rois_temp, roi_size_temp = np.unique(roi_labels, return_counts=True)

    # removing the background label
    index_bkgnd = np.argwhere(uniq_rois_temp==parcellate.null_roi_name)[0]
    uniq_rois = np.delete(uniq_rois_temp, index_bkgnd)
    roi_size  = np.delete(roi_size_temp, index_bkgnd)

    num_nodes = len(uniq_rois)

    return uniq_rois, roi_size, num_nodes


def __save(weight_vec, out_dir, subject, str_suffix = None):
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
            print('Saved the features to \n{}'.format(out_weights_path))
        except:
            print('unable to save extracted features to {}'.format(out_weights_path))
            traceback.print_exc()

    return


def __stamp_experiment(base_feature, atlas, smoothing_param, node_size, weight_method):
    "Constructs a string to uniquely identify a given feature extraction method."

    # expt_id = 'feature_{}_atlas_{}_smoothing_{}_size_{}'.format(base_feature, atlas, smoothing_param, node_size)
    expt_id = '{}_{}_smoothing{}_size{}_edgeweight_{}'.format(base_feature, atlas, smoothing_param, node_size, weight_method)

    return expt_id


def __parameter_check(base_feature, in_dir, atlas, smoothing_param, node_size, out_dir, return_results):
    """"""

    if base_feature not in __base_feature_list:
        raise NotImplementedError('Choice {} is invalid or not implemented'.format(base_feature))

    if atlas.upper() not in parcellate.atlas_list:
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


def __read_features_and_groups(features_path, groups_path):
    "Reader for data and groups"

    try:
        features = np.loadtxt(features_path)
        groups = np.loadtxt(groups_path)
    except:
        raise IOError('error reading the specified features and/or groups.')

    if len(features) != len(groups):
        raise ValueError("lengths of features and groups do not match!")

    return features, groups


def cli_run():
    "command line interface!"

    subject_ids_path, input_dir, base_feature, weight_method, \
        atlas, out_dir, node_size, smoothing_param = __parse_args()

    # when run from CLI, results will not be received
    # so no point in wasting memory maintaining a very big array
    return_results = False

    extract(subject_ids_path, input_dir, base_feature, weight_method,
            atlas, smoothing_param, node_size, out_dir, return_results)

    return


def __parse_args():
    """Parser/validator for the cmd line args."""

    parser = argparse.ArgumentParser(prog="graynet")

    parser.add_argument("-s", "--subject_ids_path", action="store", dest="subject_ids_path",
                        required=True,
                        help="Abs. path to file containing features for a given subject")

    parser.add_argument("-i", "--input_dir", action="store", dest="input_dir",
                        required=True,
                        help="path to a folder containing input data e.g. Freesurfer SUBJECTS_DIR.")

    # TODO let users specify multiple features comma separated
    parser.add_argument("-f", "--feature", action="store", dest="feature",
                        default=__default_feature, required=False,
                        help="Atlas to use to define nodes/ROIs. Default: {}".format(__default_feature))

    # TODO let users specify multiple weight methods comma separated
    parser.add_argument("-w", "--weight_method", action="store", dest="weight_method",
                        default=__default_weight_method, required=False,
                        nargs = '*', choices = __accepted_weight_list,
                        help="Method used to estimate the weight of the edge between the pair of nodes. Default : {}".format(
                            __default_weight_method))

    parser.add_argument("-a", "--atlas", action="store", dest="atlas",
                        default=__default_atlas, required=False,
                        help="Atlas to use to define nodes/ROIs. Default: {}".format(__default_atlas))

    parser.add_argument("-o", "--out_dir", action="store", dest="out_dir",
                        default=None, required=False,
                        help="Where to save the extracted features. ")

    parser.add_argument("-n", "--node_size", action="store", dest="node_size",
                        default=__default_node_size, required=False,
                        help="Parameter defining the size of individual node for the atlas parcellation. Default : {}".format(__default_node_size))

    parser.add_argument("-p", "--smoothing_param", action="store", dest="smoothing_param",
                        default=__default_smoothing_param, required=False,
                        help="Smoothing parameter for feature. "
                             "Default: FWHM of {} for Freesurfer thickness".format(__default_smoothing_param))

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

    return subject_ids_path, input_dir, params.feature, params.weight_method, \
           params.atlas, out_dir, params.node_size, params.smoothing_param


if __name__ == '__main__':
    cli_run()
