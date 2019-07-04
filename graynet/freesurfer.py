from __future__ import print_function

__all__ = ['import_features', 'get_data']

from graynet import config_graynet as cfg

from os.path import join as pjoin, exists as pexists
import collections
import nibabel
import warnings
import traceback
import numpy as np

_base_feature_list = ['thickness', 'curv', 'sulc', 'area',
                      'area.pial', 'jacobian_white']

def import_features(fs_dir,
                    subject_list,
                    base_feature= 'freesurfer_thickness',
                    fwhm=10,
                    atlas='fsaverage'):
    "Ensure subjects are provided and their data exist."

    if isinstance(subject_list, collections.Iterable):
        if len(subject_list) < 1:
            raise ValueError('Empty subject list.')
        subjects_list = subject_list
    elif isinstance(subject_list, str):
        if not pexists(subject_list):
            raise IOError('path to subject list does not exist: {}'
                          ''.format(subject_list))
        subjects_list = np.atleast_1d(np.genfromtxt(subject_list, dtype=str).astype(str))
    else:
        raise ValueError('Invalid value provided for subject list. \n '
                         'Must be a list of paths, or '
                         'path to file containing list of paths, one for each subject.')

    features= dict()
    for subject_id in subjects_list:
        try:
            print('Reading {} for {} ... '.format(base_feature, subject_id), end='')
            features[subject_id] = get_data(fs_dir, subject_id, base_feature, fwhm, atlas)
            print(' Done.')
        except:
            traceback.print_exc()
            raise ValueError('{} data for {} could not be read!'
                             ''.format(base_feature, subject_id))

    return features


def get_data(fs_dir, subject_id, base_feature, fwhm=10, atlas='fsaverage'):
    "Reads the specified features from both hemispheres for a given subject."

    feat_name = base_feature.lower()
    if feat_name in cfg.features_freesurfer:
        bare_name_feature = feat_name.replace('freesurfer_','')
        left  = __read_morph_feature(
                path_to_vertex_data(fs_dir, subject_id, hemi='lh',
                                    feature=bare_name_feature, atlas=atlas, fwhm=fwhm))
        right = __read_morph_feature(
                path_to_vertex_data(fs_dir, subject_id, hemi='rh',
                                    feature=bare_name_feature, atlas=atlas, fwhm=fwhm))
        whole = np.hstack((left, right))
    else:
        raise ValueError('Invalid choice for freesurfer data.'
                         ' Valid choices: {}'.format(cfg.features_freesurfer))

    return whole


def __all_data_exists(fs_dir, subject_id, base_feature, fwhm=10, atlas='fsaverage'):
    "Ensures all data exists for a given subject"

    if base_feature.lower() in _base_feature_list:
        data_exists = True
        for hemi in ['lh', 'rh']:
            if not pexists(path_to_vertex_data(fs_dir, subject_id,
                                               feature=base_feature,
                                               hemi=hemi,
                                               atlas=atlas, fwhm=fwhm)):
                return False
    else:
        raise ValueError('Invalid choice for freesurfer data. '
                         'Valid choices: {}'.format(_base_feature_list))

    return data_exists


def path_to_vertex_data(fsd, sid, hemi='lh', fwhm=10,
                        atlas='fsaverage', feature ='thickness'):
    "Returning the path to surface features. Using a smoothed version"

    return pjoin(fsd, sid, 'surf',
                 '{}.{}.fwhm{}.{}.mgh'.format(hemi, feature, fwhm, atlas))


def __read_morph_feature(tpath):
    "Assumes mgh format: lh.thickness.fwhm10.fsaverage.mgh"
    vec = nibabel.load(tpath).get_data() # typically of shape: (163842, 1, 1)

    return np.squeeze(vec) # becomes (163842, )


