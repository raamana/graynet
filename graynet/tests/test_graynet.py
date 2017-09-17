import os
import shlex
import sys
from os.path import join as pjoin, exists as pexists, abspath
from sys import version_info
import numpy as np
import scipy.stats
sys.dont_write_bytecode = True

from pytest import raises, warns, set_trace

if version_info.major==2 and version_info.minor==7:
    import graynet
    from graynet import cli_run as CLI
elif version_info.major > 2:
    from graynet import run_workflow
    from graynet import cli_run as CLI
else:
    raise NotImplementedError('hiwenet supports only 2.7.13 or 3+. Upgrate to Python 3+ is recommended.')


test_dir = os.path.dirname(os.path.realpath(__file__))
base_dir = os.path.realpath( pjoin(test_dir, '..', '..', 'example_data') )

subject_id_list = ['subject12345', ]

out_dir = pjoin(base_dir, 'graynet')
if not pexists(out_dir):
    os.mkdir(out_dir)

fs_dir = pjoin(base_dir, 'freesurfer')
base_feature = 'freesurfer_thickness'
atlas = 'FSAVERAGE' # 'GLASSER2016' #
fwhm = 10

num_roi_atlas = {'FSAVERAGE':68, 'GLASSER2016':360}
num_roi_wholebrain = num_roi_atlas[atlas]
num_links = num_roi_wholebrain*(num_roi_wholebrain-1)/2

weight_methods = [ 'manhattan', 'minowski' ]

cur_dir = os.path.dirname(abspath(__file__))
example_dir = abspath(pjoin(cur_dir, '..', '..', 'example_data', 'freesurfer'))
sub_list = pjoin(example_dir, 'subject_list.txt')
out_dir = pjoin(example_dir, 'test_outputs')

def test_run_no_IO():
    edge_weights_all = graynet.extract(subject_id_list, fs_dir, base_feature, weight_methods, atlas, fwhm,
                                       out_dir=None, return_results=True)
    num_combinations = len(list(edge_weights_all))

    if num_combinations != len(subject_id_list)*len(weight_methods):
        raise ValueError('invalid results : # subjects')

    for wm in weight_methods:
        for sub in subject_id_list:
            if edge_weights_all[(wm, sub)].size != num_links:
                raise ValueError('invalid results : # links')

def test_run_roi_stats():
    "Tests whether roi stats can be computed (not their accuracy) and the return values match in size."

    summary_methods = ['median', 'mode', 'mean', 'std', 'gmean', 'hmean', 'variation', 'entropy', 'skew', 'kurtosis']
    trimmed_mean = lambda array: scipy.stats.trim_mean(array, proportiontocut=0.05)
    third_kstat = lambda array: scipy.stats.kstat(array, n=3)
    summary_methods.extend([trimmed_mean, third_kstat])
    # checking support for nan-handling callables
    summary_methods.extend([np.nanmedian, np.nanmean])

    for summary_method in summary_methods:
        roi_medians = graynet.roiwise_stats_indiv(subject_id_list, fs_dir, base_feature,
                                                  summary_method, atlas, fwhm, out_dir=None, return_results=True)
        for sub in subject_id_list:
            if roi_medians[sub].size != num_roi_wholebrain:
                raise ValueError('invalid summary stats - #nodes do not match.')


def test_CLI_weight():
    " ensures the CLI works. "

    sys.argv = shlex.split('graynet -s {} -i {} -w manhattan -o {} -a {}'.format(sub_list, example_dir, out_dir, atlas))

    CLI()

def test_CLI_stats():
    " ensures the CLI works. "

    sys.argv = shlex.split('graynet -s {} -i {} -r median gmean -o {} -a {}'.format(sub_list, example_dir, out_dir, atlas))

    CLI()

def test_CLI_only_weight_or_stats():
    " ensures the CLI works. "

    with warns(UserWarning):
        sys.argv = shlex.split('graynet -s {} -i {} -w cosine -r median gmean -o {} -a {}'.format(sub_list, example_dir, out_dir, atlas))

    CLI()

test_CLI_stats()