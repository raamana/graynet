
import sys
import os
import nibabel
import warnings
import networkx as nx
import numpy as np

_base_feature_list = ['thickness', 'gmdensity']

def extract(fs_dir, subject_id_list, base_feature, atlas ='fsaverage', fwhm = 10, size = 100):
    """Extracts weighted networks from gray matters features based on Freesurfer processing."""

    __parameter_check(fs_dir, base_feature, atlas, fwhm, size)

    num_nodes = 10
    ewmat = np.zeros([num_nodes, num_nodes])

    return ewmat


def __parameter_check(fs_dir, base_feature, atlas, fwhm, size):
    """"""

    if base_feature not in _base_feature_list:
        raise NotImplementedError('Invalid choice, or network extraction for {} is not implemented'.format(base_feature))


if __name__ == '__main__':
    extract()