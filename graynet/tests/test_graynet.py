import os
import shlex
import sys
from os.path import join as pjoin, exists as pexists, abspath, dirname, realpath
from sys import version_info
import numpy as np
import scipy.stats

sys.dont_write_bytecode = True

from pytest import raises, warns, set_trace

if __name__ == '__main__' and __package__ is None:
    parent_dir = dirname(dirname(abspath(__file__)))
    pkg_dir = dirname(parent_dir)
    sys.path.append(parent_dir)
    sys.path.append(pkg_dir)

if version_info.major > 2:
    from graynet import run_workflow as graynet
    from graynet import cli_run as CLI
    from graynet import multi_edge
    from graynet import config_graynet as cfg
else:
    raise NotImplementedError('graynet requires Python 3+.')

test_dir = dirname(os.path.realpath(__file__))
base_dir = realpath(pjoin(test_dir, '..', '..', 'example_data'))

subject_id_list = ['subject12345', ]

out_dir = pjoin(base_dir, 'graynet')
if not pexists(out_dir):
    os.mkdir(out_dir)

fs_dir = pjoin(base_dir, 'freesurfer')
base_feature = 'freesurfer_thickness'
atlas = 'FSAVERAGE'  # 'GLASSER2016' #
fwhm = 10

num_roi_atlas = {'FSAVERAGE': 68, 'GLASSER2016': 360}
num_roi_wholebrain = num_roi_atlas[atlas]
num_links = num_roi_wholebrain * (num_roi_wholebrain - 1) / 2

weight_methods = ['manhattan', ]

cur_dir = os.path.dirname(abspath(__file__))
example_dir = abspath(pjoin(cur_dir, '..', '..', 'example_data', 'freesurfer'))
sub_list = pjoin(example_dir, 'subject_list.txt')
out_dir = pjoin(example_dir, 'test_outputs')
if not pexists(out_dir):
    os.mkdir(out_dir)

dimensionality = 1000
num_groups = 5

cur_dir = os.path.dirname(abspath(__file__))


def test_multi_edge():
    edge_weights_all = multi_edge.extract_multiedge(subject_id_list,
                                                    input_dir=fs_dir,
                                                    base_feature_list=cfg.default_features_multi_edge,
                                                    edge_range_dict=cfg.edge_range_predefined,
                                                    weight_method_list=weight_methods,
                                                    atlas=atlas,
                                                    smoothing_param=fwhm,
                                                    out_dir=out_dir,
                                                    return_results=True,
                                                    num_procs=1)

    num_combinations = len(list(edge_weights_all))
    expected_num_comb = len(subject_id_list) * len(weight_methods)*len(cfg.default_features_multi_edge)
    if num_combinations != expected_num_comb:
        raise ValueError('invalid results : # subjects')

    for wm in weight_methods:
        for sub in subject_id_list:
            for feat in cfg.default_features_multi_edge:
                if edge_weights_all[(wm, feat, sub)].size != num_links:
                    raise ValueError('invalid results : # links')

    print('')

def test_multi_edge_CLI():

    sys.argv = shlex.split('graynet -s {} -i {} '
                           ' -f freesurfer_thickness freesurfer_curv'
                           ' --do_multi_edge --multi_edge_range 0.0 5.0 -0.3 +0.3 '
                           ' -w manhattan -o {} -a {}'.format(sub_list, example_dir, out_dir, atlas))

    CLI()


def test_run_no_IO():
    edge_weights_all = graynet.extract(subject_id_list,
                                       fs_dir,
                                       base_feature=base_feature,
                                       weight_method_list= weight_methods,
                                       atlas=atlas,
                                       smoothing_param=fwhm,
                                       out_dir=out_dir,
                                       return_results=True,
                                       num_procs=4)
    num_combinations = len(list(edge_weights_all))

    if num_combinations != len(subject_id_list) * len(weight_methods):
        raise ValueError('invalid results : # subjects')

    for wm in weight_methods:
        for sub in subject_id_list:
            if edge_weights_all[(wm, sub)].size != num_links:
                raise ValueError('invalid results : # links')


def test_run_roi_stats_via_API():
    "Tests whether roi stats can be computed (not their accuracy) and the return values match in size."

    summary_methods = ['median', 'mean', 'std', 'variation', 'entropy', 'skew', 'kurtosis']
    # 'mode' returns more than one value; 'gmean' requires only positive values,
    # 'hmean' can not always be computed
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


def test_run_roi_stats_via_CLI():
    " ensures the CLI works. "

    sys.argv = shlex.split(
        'graynet -s {} -i {} -r median gmean -o {} -a {}'.format(sub_list, example_dir, out_dir, atlas))

    CLI()


def test_CLI_only_weight_or_stats():
    " ensures the CLI works. "

    with raises(SystemExit):
        sys.argv = shlex.split(
            'graynet -s {} -i {} -w cosine -r median gmean -o {} -a {}'.format(sub_list, example_dir, out_dir, atlas))
        CLI()


def test_empty_subject_list():
    with raises(ValueError):
        ew = graynet.extract([], fs_dir)


def test_invalid_edge_range():
    with raises(ValueError):
        ew = graynet.extract(subject_id_list, fs_dir, edge_range=-1)

    with raises(ValueError):
        ew = graynet.extract(subject_id_list, fs_dir, edge_range=[])

    with raises(ValueError):
        ew = graynet.extract(subject_id_list, fs_dir, edge_range=[1, ])

    with raises(ValueError):
        ew = graynet.extract(subject_id_list, fs_dir, edge_range=[1, 2, 3])

    with raises(ValueError):
        ew = graynet.extract(subject_id_list, fs_dir, edge_range=(1, np.NaN))

    with raises(ValueError):
        ew = graynet.extract(subject_id_list, fs_dir, edge_range=(2, 1))


def test_invalid_nbins():
    with raises(ValueError):
        ew = graynet.extract(subject_id_list, fs_dir, num_bins=np.NaN)

    with raises(ValueError):
        ew = graynet.extract(subject_id_list, fs_dir, num_bins=np.Inf)

    with raises(ValueError):
        ew = graynet.extract(subject_id_list, fs_dir, num_bins=2)


#
test_multi_edge_CLI()

# test_empty_subject_list()
# test_run_no_IO()
# test_run_roi_stats_via_API()
# test_run_roi_stats_via_CLI()
# test_CLI_only_weight_or_stats()
