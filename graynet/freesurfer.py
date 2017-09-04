
from os.path import join as pjoin, exists as pexists
import collections
import nibabel
import warnings
import numpy as np


def import_features(fs_dir, subject_list, base_feature= 'thickness', fwhm=10, atlas='fsaverage'):
    "Ensure subjects are provided and their data exist."

    if isinstance(subject_list, collections.Iterable):
        if len(subject_list) < 1:
            raise ValueError('Empty subject list.')
        subjects_list = subject_list
    elif isinstance(subject_list, str):
        if not pexists(subject_list):
            raise IOError('path to subject list does not exist: {}'.format(subject_list))
        subjects_list = np.loadtxt(subject_list, dtype=str)
    else:
        raise ValueError('Invalid value provided for subject list. \n '
                         'Must be a list of paths, or path to file containing list of paths, one for each subject.')

    features= dict()
    for subject_id in subjects_list:
        try:
            features[subject_id] = __get_data(fs_dir, subject_id, base_feature, fwhm, atlas)
        except:
            raise ValueError('{} data for {} could not be read!'.format(base_feature, subject_id))

    return features


def __get_data(fs_dir, subject_id, base_feature, fwhm=10, atlas='fsaverage'):
    "Ensures all data exists for a given subject"

    _base_feature_list = ['freesurfer_thickness', 'freesurfer_curv', 'freesurfer_sulc']
    feat_name = base_feature.lower()
    if feat_name in _base_feature_list:
        bare_name_feature = feat_name.replace('freesurfer_','')
        left  = __read_morph_feature(_surf_data_path(fs_dir, subject_id, hemi='lh', feature=bare_name_feature))
        right = __read_morph_feature(_surf_data_path(fs_dir, subject_id, hemi='rh', feature=bare_name_feature))
        whole = np.hstack((left, right))
    else:
        raise ValueError('Invalid choice for freesurfer data. Valid choices: {}'.format(_base_feature_list))

    return whole


def __all_data_exists(fs_dir, subject_id, base_feature):
    "Ensures all data exists for a given subject"

    _base_feature_list = ['thickness', 'curv', 'sulc']
    if base_feature.lower() in _base_feature_list:
        data_exists = pexists(_surf_data_path(fs_dir, subject_id, 'lh')) and \
                      pexists(_surf_data_path(fs_dir, subject_id, 'rh'))
    else:
        raise ValueError('Invalid choice for freesurfer data. Valid choices: {}'.format(_base_feature_list))

    return data_exists


def _surf_data_path(fsd, sid, hemi='lh', fwhm=10, atlas='fsaverage', feature = 'thickness'):
    "Returning the path to surface features. Using a smoothed version"

    return pjoin(fsd, sid, 'surf', '{}.{}.fwhm{}.{}.mgh'.format(hemi, feature, fwhm, atlas))


def __read_morph_feature(tpath):
    "Assumes mgh format: lh.thickness.fwhm10.fsaverage.mgh"
    vec = nibabel.load(tpath).get_data() # typically of shape: (163842, 1, 1)

    return np.squeeze(vec) # becomes (163842, )


def __read_data(fs_dir, subject_list, base_feature):
    "Returns the location of the source of subject-wise features: /path/subject/surf/?h.thickness or nifti image"

    def read_gmdensity(gmpath):
        return nibabel.load(gmpath)

    reader = {'gmdensity': read_gmdensity, 'thickness': __read_morph_feature}

    features = dict()
    for subj_info in subject_list:
        subj_id, subj_data_path = subj_info.split(',')
        try:
            features[subj_id] = reader[base_feature](subj_data_path)
        except:
            warnings.warn('data for {} could not be read from:\n{}'.format(subj_id, subj_data_path))
            raise

    return features
