import os
from os.path import join as pjoin, exists as pexists

from graynet import extract

# base_dir = '/u1/work/hpc3194/PPMI/processed'
# subject_id_list = ['4139_bl_PPMI', '4118_bl_PPMI', '4116_bl_PPMI', '4105_bl_PPMI']

test_dir = os.path.dirname(os.path.realpath(__file__))
base_dir = os.path.realpath( pjoin(test_dir, '..', 'example_data') )

subject_id_list = ['4113_bl_PPMI', ]

out_dir = pjoin(base_dir, 'graynet')
if not pexists(out_dir):
    os.mkdir(out_dir)

fs_dir = pjoin(base_dir, 'freesurfer')
base_feature = 'thickness'
atlas ='GLASSER2016'
fwhm = 10

extract(subject_id_list, fs_dir, base_feature, atlas, fwhm)

