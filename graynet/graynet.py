import collections
import warnings
import os
from os.path import join as pjoin, exists as pexists

import hiwenet
import nibabel
import numpy as np
import traceback

import freesurfer
import parcellate

np.seterr(divide='ignore', invalid='ignore')

__base_feature_list = ['thickness', 'gmdensity']

__default_weight_method = 'minowski' # 'manhattan'
__default_num_bins = 25
__default_trim_percentile = 5


def extract(subject_id_list, fs_dir,
            base_feature='thickness',
            weight_method= __default_weight_method,
            atlas ='GLASSER2016',
            fwhm = 10, size = None,
            out_dir=None):
    """Extracts weighted networks from gray matters features based on Freesurfer processing. Subject_id_list must be a file or list containing one id,path for each subject. """

    __parameter_check(base_feature, fs_dir, atlas, fwhm, size)
    subject_id_list = __subject_check(subject_id_list)
    num_subjects = len(subject_id_list)

    roi_labels, ctx_annot = parcellate.freesurfer_roi_labels(atlas)
    uniq_rois, roi_size, num_nodes = __roi_info(roi_labels)

    features = freesurfer.import_features(fs_dir, subject_id_list, base_feature)

    edge_weights_all = np.zeros([num_subjects, num_nodes*(num_nodes-1)/2])
    for ss, subject in enumerate(subject_id_list):

        print('Processing {} : {} / {}'.format(subject, ss+1, num_subjects))

        try:
            data, rois = __remove_background_roi(features[subject], roi_labels, parcellate.null_roi_index)
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


def __read_data(subject_list, base_feature):
    "Returns the location of the source of subject-wise features: /path/subject/surf/?h.thickness or nifti image"

    def read_thickness(tpath):
        return nibabel.freesurfer.io.read_morph_data(tpath)

    def read_gmdensity(gmpath):
        return nibabel.load(gmpath)

    reader = {'gmdensity': read_gmdensity, 'thickness': read_thickness}

    features = dict()
    for subj_info in subject_list:
        subj_id, subj_data_path = subj_info.split(',')
        try:
            features[subj_id] = reader[base_feature](subj_data_path)
        except:
            warnings.warn('data for {} could not be read from:\n{}'.format(subj_id, subj_data_path))
            raise

    return features


def __remove_background_roi(data,labels, ignore_label):
    "Returns everything but specified label"

    mask = labels != ignore_label

    return data[mask], labels[mask]


def __roi_info(roi_labels):
    "Unique ROIs in a given atlas parcellation, count and size. Excludes the background"

    uniq_rois_temp, roi_size = np.unique(roi_labels, return_counts=True)

    # removing the background label
    background_labels = [ parcellate.null_roi_index, ]
    uniq_rois = np.delete(uniq_rois_temp, background_labels)

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


def __parameter_check(base_feature, fs_dir, atlas, fwhm, size):
    """"""

    if base_feature not in __base_feature_list:
        raise NotImplementedError('Choice {} is invalid or not implemented'.format(base_feature))

    if atlas.upper() not in parcellate.atlas_list:
        raise ValueError('Invalid atlas choice. Use one of {}'.format(parcellate.atlas_list))

    if not pexists(fs_dir):
        raise IOError('Freesurfer directory at {} does not exist.'.format(fs_dir))

    # no checks on subdivison size yet, as its not implemented

    return


def __cli_run():
    "command line interface!"

    raise NotImplementedError('command line interface not implemented yet.')


if __name__ == '__main__':
    __cli_run()