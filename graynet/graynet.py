import collections
import os
import traceback
import warnings
from os.path import join as pjoin, exists as pexists

import hiwenet
import nibabel
import numpy as np

import freesurfer
import parcellate

np.seterr(divide='ignore', invalid='ignore')

__features_freesurfer = ['freesurfer_thickness', ]
__features_fsl = ['gmdensity', ]

__base_feature_list = __features_freesurfer + __features_fsl

__default_weight_method = 'minowski' # 'manhattan'
__default_num_bins = 25
__default_trim_percentile = 5


def extract(subject_id_list, input_dir, base_feature = 'thickness', weight_method = __default_weight_method,
            atlas ='GLASSER2016', smoothing_param = 10, node_size = None, out_dir=None):
    """
    Extracts weighted networks from gray matters features based on Freesurfer processing.
    Subject_id_list must be a file or a list containing one id,path for each subject.

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
    weight_method : str
        Name of covariance metric to use to compute weights. Currently those supported by hiwenet, which can be one of:
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

    Returns
    -------
    edge_weights_all : ndarray
        Numpy array of size n x p, where n = number of subjects and p = k*(k-1)/2,
        where k = number of nodes in the atlas parcellation.

    """

    __parameter_check(base_feature, input_dir, atlas, smoothing_param, node_size)
    subject_id_list = __subject_check(subject_id_list)
    num_subjects = len(subject_id_list)

    roi_labels, ctx_annot = parcellate.freesurfer_roi_labels(atlas)
    uniq_rois, roi_size, num_nodes = __roi_info(roi_labels)

    features = import_features(input_dir, subject_id_list, base_feature)

    edge_weights_all = np.zeros([num_subjects, num_nodes*(num_nodes-1)/2])
    for ss, subject in enumerate(subject_id_list):

        print('Processing {} : {} / {}'.format(subject, ss+1, num_subjects))

        try:
            data, rois = __remove_background_roi(features[subject], roi_labels, parcellate.null_roi_name)
            edge_weights = hiwenet.extract(data, rois, weight_method)

            weight_vec = __get_triu_handle_inf_nan(edge_weights)

            __save(weight_vec, out_dir, subject)
            edge_weights_all[ss, :] = weight_vec

        except (RuntimeError, RuntimeWarning) as runexc:
            print(runexc)
            pass
        except:
            print('Unable to extract covariance features for {}'.format(subject))
            traceback.print_exc()

    return edge_weights_all


def import_features(input_dir, subject_id_list, base_feature):
    "Wrapper to support input data of multiple types and multiple packages."

    base_feature = base_feature.lower()
    if base_feature in __features_freesurfer:
        features = freesurfer.import_features(input_dir, subject_id_list, base_feature)
    elif base_feature in __features_fsl:
        features = fsl_import(input_dir, subject_id_list, base_feature)
    else:
        raise ValueError('Invalid choice. Choose one of \n {}'.format(__base_feature_list))

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
        warnings.warn(' {} non-finite values are found.'.format(num_nonfinite))

    return upper_tri_vec


def __subject_check(subjects_info):
    "Ensure subjects are provided and their data exist."

    if isinstance(subjects_info, collections.Iterable):
        if len(subjects_info) < 1:
            raise ValueError('Empty subject list.')
        subjects_list = subjects_info
    elif isinstance(subjects_info, str):
        if not pexists(subjects_info):
            raise IOError('path to subject list does not exist: {}'.format(subjects_info))
        subjects_list = np.loadtxt(subjects_info, dtype=str)
    else:
        raise ValueError('Invalid value provided for subject list. \n '
                         'Must be a list of paths, or path to file containing list of paths, one for each subject.')

    return subjects_list


def __remove_background_roi(data,labels, ignore_label):
    "Returns everything but specified label"

    mask = labels != ignore_label

    return data[mask], labels[mask]


def __roi_info(roi_labels):
    "Unique ROIs in a given atlas parcellation, count and size. Excludes the background"

    uniq_rois_temp, roi_size = np.unique(roi_labels, return_counts=True)

    # removing the background label
    index_bkgnd = np.argwhere(uniq_rois_temp==parcellate.null_roi_name)[0]
    uniq_rois = np.delete(uniq_rois_temp, index_bkgnd)

    num_nodes = len(uniq_rois)

    return uniq_rois, roi_size, num_nodes


def __save(weight_vec, out_dir, subject):
    "Saves the features to disk."

    if out_dir is not None:
        # get outpath returned from hiwenet, based on dist name and all other parameters
        # choose out_dir name  based on dist name and all other parameters
        out_subject_dir = pjoin(out_dir, subject)
        if not pexists(out_subject_dir):
            os.mkdir(out_subject_dir)
        out_weights_path = pjoin(out_subject_dir, 'graynet.csv')

        try:
            np.savetxt(out_weights_path, weight_vec, fmt='%.5f')
            print('Saved the features to \n{}'.format(out_weights_path))
        except:
            print('unable to save extracted features to {}'.format(out_weights_path))
            traceback.print_exc()

    return


def __parameter_check(base_feature, in_dir, atlas, smoothing_param, node_size):
    """"""

    if base_feature not in __base_feature_list:
        raise NotImplementedError('Choice {} is invalid or not implemented'.format(base_feature))

    if atlas.upper() not in parcellate.atlas_list:
        raise ValueError('Invalid atlas choice. Use one of {}'.format(parcellate.atlas_list))

    if not pexists(in_dir):
        raise IOError('Input directory at {} does not exist.'.format(in_dir))

    # no checks on subdivison size yet, as its not implemented

    return


def __cli_run():
    "command line interface!"

    raise NotImplementedError('command line interface not implemented yet.')


if __name__ == '__main__':
    __cli_run()