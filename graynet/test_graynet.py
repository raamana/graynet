import os
from os.path import join as pjoin, exists as pexists
from graynet import extract

test_dir = os.path.dirname(os.path.realpath(__file__))
base_dir = os.path.realpath( pjoin(test_dir, '..', 'example_data') )

subject_id_list = ['subject12345', ]

out_dir = pjoin(base_dir, 'graynet')
if not pexists(out_dir):
    os.mkdir(out_dir)

fs_dir = pjoin(base_dir, 'freesurfer')
base_feature = 'thickness'
atlas ='GLASSER2016'
fwhm = 10
num_roi_wholebrain = 360
num_links = num_roi_wholebrain*(num_roi_wholebrain-1)/2

weight_method = 'manhattan' # 'minowski' # 'manhattan'

def test_run_no_IO():
    edge_weights_all = extract(subject_id_list, fs_dir, base_feature, weight_method, atlas, fwhm)

    ew_shape = edge_weights_all.shape
    assert ew_shape[0] == len(subject_id_list)
    assert ew_shape[1] == num_links

