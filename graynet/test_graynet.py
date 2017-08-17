

import sys
import os
from os.path import join as pjoin, exists as pexists

import graynet

base_dir = '/u1/work/hpc3194/PPMI/processed'

subject_id_list = ['4139_bl_PPMI', '4118_bl_PPMI', '4116_bl_PPMI', '4105_bl_PPMI']
out_dir = pjoin(base_dir, 'graynet')
if not pexists(out_dir):
    os.mkdir(out_dir)

fs_dir = pjoin(base_dir, 'freesurfer')
base_feature = 'thickness'
atlas ='GLASSER2016'
fwhm = 10

graynet.extract(subject_id_list, out_dir, fs_dir, base_feature, atlas, fwhm)

