
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

_base_feature_list = ['thickness', 'gmdensity']

def extract(subject_id_list, base_feature, atlas ='GLASSER2016', fwhm = 10, size = None):
    """Extracts weighted networks from gray matters features based on Freesurfer processing."""

    __parameter_check(base_feature, atlas, fwhm, size)

    __subject_check(subject_id_list, base_feature)

    num_nodes = 10
    ewmat = np.zeros([num_nodes, num_nodes])

    return ewmat


def __subject_check(subjects_info, base_feature):
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

    nonexisting_subjects = [ data_path for data_path in subjects_list if not pexists(data_path)]
    if len(nonexisting_subjects) > 0:
        raise ValueError('Following {} subjects do not exist:\n {}'.format(
            len(nonexisting_subjects), '\n'.join(nonexisting_subjects)))

    return


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


def __parameter_check(base_feature, atlas, fwhm, size):
    """"""

    if base_feature not in _base_feature_list:
        raise NotImplementedError('Choice {} is invalid or not implemented'.format(base_feature))

    if atlas.upper() not in parcellate.atlas_list:
        raise ValueError('Invalid atlas choice. Use one of {}'.format(parcellate.atlas_list))

    # no checks on subdivison size yet, as its not implemented

    return


def cli_run():
    "command line interface!"

    raise NotImplementedError('command line interface not implemented yet.')


if __name__ == '__main__':
    cli_run()