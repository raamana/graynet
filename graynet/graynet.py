
import sys
import os
from os.path import join as pjoin, exists as pexists
import collections
import nibabel
import warnings
import networkx as nx
import numpy as np
import hiwenet
from . import parcellate
from . import freesurfer

_base_feature_list = ['thickness', 'gmdensity']

default_weight_method = 'manhattan'
default_num_bins = 25
default_trim_percentile = 5


def extract(subject_id_list, out_dir, fs_dir,
            base_feature='thickness', atlas ='GLASSER2016',
            fwhm = 10, size = None):
    """Extracts weighted networks from gray matters features based on Freesurfer processing. Subject_id_list must be a file or list containing one id,path for each subject. """

    __parameter_check(base_feature, fs_dir, atlas, fwhm, size)
    subject_id_list = __subject_check(subject_id_list)
    num_subjects = len(subject_id_list)

    ctx_parc = parcellate.get_atlas_annot(atlas)
    num_nodes = __num_nodes(ctx_parc)

    features = freesurfer.import_features(fs_dir, subject_id_list, base_feature)

    ewmat = np.zeros([num_subjects, num_nodes])
    for subject in subject_id_list:
        out_weights_path= pjoin(out_dir,subject, 'graynet.csv')
        try:
            edge_weights = hiwenet.extract(features[subject], ctx_parc,
                                           weight_method= default_weight_method,
                                           num_bins=default_num_bins, trim_outliers=True,
                                           trim_percentile=default_trim_percentile)
            weights_array = edge_weights[ np.triu_indices_from(edge_weights, 1) ]
        except:
            raise ValueError('Unable to extract covariance features for {}'.format(subject))

        try:
            np.save(weights_array, out_weights_path)
        except:
            raise IOError('unable to save extracted features to {}'.format(out_weights_path))


    return ewmat


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


def __num_nodes(atlas_annot):
    "Number of unique ROIs in a given atlas parcellation"


def __parameter_check(base_feature, fs_dir, atlas, fwhm, size):
    """"""

    if base_feature not in _base_feature_list:
        raise NotImplementedError('Choice {} is invalid or not implemented'.format(base_feature))

    if atlas.upper() not in parcellate.atlas_list:
        raise ValueError('Invalid atlas choice. Use one of {}'.format(parcellate.atlas_list))

    if not pexists(fs_dir):
        raise IOError('Freesurfer directory at {} does not exist.'.format(fs_dir))

    # no checks on subdivison size yet, as its not implemented

    return


def cli_run():
    "command line interface!"

    raise NotImplementedError('command line interface not implemented yet.')


if __name__ == '__main__':
    cli_run()