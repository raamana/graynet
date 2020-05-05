import os
import shlex
import sys
from os.path import abspath, dirname, exists as pexists, join as pjoin, realpath
from sys import version_info

import numpy as np

sys.dont_write_bytecode = True

import traceback
from pytest import raises
from hypothesis import given, strategies
from hypothesis import settings as hyp_settings

if __name__ == '__main__' and __package__ is None:
    parent_dir = dirname(dirname(abspath(__file__)))
    pkg_dir = dirname(parent_dir)
    sys.path.append(parent_dir)
    sys.path.append(pkg_dir)

if version_info.major > 2:
    from graynet import config_graynet as cfg
    from graynet.run_workflow import cli_run as CLI
    from graynet import run_workflow as graynet
    from graynet.run_workflow import extract
    from graynet.multi_edge import extract_multiedge
else:
    raise NotImplementedError('graynet requires Python 3+.')

test_dir = dirname(os.path.realpath(__file__))
base_dir = realpath(pjoin(test_dir, '..', '..', 'example_data'))

out_dir = pjoin(base_dir, 'graynet')
if not pexists(out_dir):
    os.mkdir(out_dir)


fs_dir = pjoin(base_dir, 'freesurfer')
subject_id_list = ['subject12345', ]

base_feature = 'freesurfer_thickness'
atlas = 'fsaverage'  # 'glasser2016' #
fwhm = 10

vbm_in_dir = pjoin(base_dir, 'volumetric_CAT12')
vbm_sub_list = ['CAM_0002_01', ]

base_feature_list = ('freesurfer_thickness',
                     'spm_cat_gmdensity')
num_base_features = len(base_feature_list)

# 'glasser2016' not tested regularly
feature_to_atlas_list = {'freesurfer_thickness': ('fsaverage', ),
                         'spm_cat_gmdensity'   : (
                         'cat_aal', 'cat_lpba40', 'cat_ibsr')}

feature_to_in_dir = {'freesurfer_thickness': fs_dir,
                     'spm_cat_gmdensity'   : vbm_in_dir}
feature_to_subject_id_list = {'freesurfer_thickness': subject_id_list,
                              'spm_cat_gmdensity'   : vbm_sub_list}

num_roi_atlas = {'fsaverage'  : 68,
                 'glasser2016': 360,
                 'cat_aal'    : 122,
                 'cat_lpba40' : 56,
                 'cat_ibsr'   : 32}
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

# TODO tests for volumetric version of multiedge to be done!
def test_multi_edge():
    edge_weights_all = extract_multiedge(subject_id_list,
                                         input_dir=fs_dir,
                                         base_feature_list=cfg.default_features_multi_edge,
                                         edge_range_dict=cfg.edge_range_predefined,
                                         weight_method_list=weight_methods,
                                         atlas=atlas,
                                         smoothing_param=fwhm,
                                         out_dir=out_dir,
                                         return_results=True,
                                         num_procs=1,
                                         overwrite_results=True)

    num_combinations = len(list(edge_weights_all))
    expected_num_comb = len(subject_id_list)*len(weight_methods)*len(cfg.default_features_multi_edge)
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
                           ' -w manhattan -o {} -a {}'
                           ''.format(sub_list, example_dir, out_dir, atlas))

    CLI()

def test_multi_edge_summary_stat_CLI():

    ss_list = ' '.join(['prod', 'median', 'max', 'min', 'gmean', 'hmean', 'std'])
    sys.argv = shlex.split('graynet -s {} -i {} '
                           ' -f freesurfer_thickness freesurfer_curv'
                           ' --do_multi_edge --multi_edge_range 0.0 5.0 -0.3 +0.3 '
                           ' -w manhattan cosine --summary_stat {} '
                           '-o {} -a {}'
                           ''.format(sub_list, example_dir, ss_list, out_dir, atlas))

    CLI()

@hyp_settings(max_examples=num_base_features, deadline=None)
@given(strategies.sampled_from(base_feature_list))
def test_run_no_IO(base_feature):

    for atlas in feature_to_atlas_list[base_feature]:
        try:
            sud_id_list = feature_to_subject_id_list[base_feature]
            edge_weights_all = graynet.extract(sud_id_list,
                                               feature_to_in_dir[base_feature],
                                               base_feature=base_feature,
                                               weight_method_list= weight_methods,
                                               atlas=atlas,
                                               smoothing_param=fwhm,
                                               out_dir=out_dir,
                                               return_results=True,
                                               num_procs=1)
            num_combinations = len(list(edge_weights_all))

            if num_combinations != len(sud_id_list) * len(weight_methods):
                raise ValueError('invalid results : # subjects')

            num_roi_wholebrain = num_roi_atlas[atlas]
            num_links = num_roi_wholebrain * (num_roi_wholebrain - 1) / 2

            for wm in weight_methods:
                for sub in sud_id_list:
                    if edge_weights_all[(wm, sub)].size != num_links:
                        raise ValueError('invalid results : # links')
        except:
            traceback.print_exc()
            raise

@hyp_settings(max_examples=num_base_features, deadline=None)
@given(strategies.sampled_from(base_feature_list))
def test_run_API_on_original_features(base_feature):

    for atlas in feature_to_atlas_list[base_feature]:
        sud_id_list = feature_to_subject_id_list[base_feature]
        edge_weights_all = extract(sud_id_list,
                                   feature_to_in_dir[base_feature],
                                   base_feature=base_feature,
                                   weight_method_list=cfg.weights_on_original_features,
                                   atlas=atlas,
                                   smoothing_param=fwhm,
                                   out_dir=out_dir,
                                   return_results=True,
                                   num_procs=1)

        num_combinations = len(list(edge_weights_all))

        if num_combinations != len(sud_id_list) * len(cfg.weights_on_original_features):
            raise ValueError('invalid results : # subjects')

        num_roi_wholebrain = num_roi_atlas[atlas]
        num_links = num_roi_wholebrain * (num_roi_wholebrain - 1) / 2

        for wm in cfg.weights_on_original_features:
            for sub in sud_id_list:
                if edge_weights_all[(wm, sub)].size != num_links:
                    raise ValueError('invalid results : # links')

@hyp_settings(max_examples=num_base_features, deadline=None)
@given(strategies.sampled_from(base_feature_list))
def test_run_roi_stats_via_API(base_feature):
    """Tests whether roi stats can be computed (not their accuracy)
    and the return values match in size."""

    summary_methods = ['median', 'mean', 'std', 'variation', 'entropy', 'skew',
                       'kurtosis']
    # 'mode' returns more than one value; 'gmean' requires only positive values,
    # 'hmean' can not always be computed
    from scipy.stats import  trim_mean, kstat
    from functools import partial
    trimmed_mean = partial(trim_mean, proportiontocut=0.05)
    third_kstat = partial(kstat, n=3)

    summary_methods.extend([trimmed_mean, third_kstat])
    # checking support for nan-handling callables
    summary_methods.extend([np.nanmedian, np.nanmean])

    sud_id_list = feature_to_subject_id_list[base_feature]
    for atlas in feature_to_atlas_list[base_feature]:
        num_roi_wholebrain = num_roi_atlas[atlas]
        for summary_method in summary_methods:
            roi_medians = graynet.roiwise_stats_indiv(sud_id_list,
                                                      feature_to_in_dir[base_feature],
                                                      base_feature=base_feature,
                                                      chosen_roi_stats=summary_method,
                                                      atlas=atlas,
                                                      smoothing_param=fwhm,
                                                      out_dir=out_dir,
                                                      return_results=True)

            for sub in sud_id_list:
                if roi_medians[sub].size != num_roi_wholebrain:
                    raise ValueError('invalid summary stats - #nodes do not match.')


def test_CLI_weight():
    " ensures the CLI works. "

    sys.argv = shlex.split('graynet -s {} -i {} -w manhattan -o {} -a {}'
                           ''.format(sub_list, example_dir, out_dir, atlas))

    CLI()


def test_run_roi_stats_via_CLI():
    " ensures the CLI works. "

    sys.argv = shlex.split('graynet -s {} -i {} -r median gmean -o {} -a {}'
                           ''.format(sub_list, example_dir, out_dir, atlas))

    CLI()


def test_CLI_only_weight_or_stats():
    " ensures the CLI works. "

    with raises(SystemExit):
        sys.argv = shlex.split(
            'graynet -s {} -i {} -w cosine -r median gmean -o {} -a {}'
            ''.format(sub_list, example_dir, out_dir, atlas))
        CLI()


def test_empty_subject_list():
    # API
    with raises(ValueError):
        ew = graynet.extract([], fs_dir)

    # in CLI, only non-Freesurfer lead to an error
    for feat in cfg.features_volumetric:  # invalid list
        with raises(ValueError):
            sys.argv = shlex.split('graynet -i {} -f {}'.format(fs_dir, feat))
            run_cli()


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


def test_atlas_parcel_subdivision():

    wm = 'manhattan'
    # much slower: zip(cfg.allowed_mvpp, cfg.mvpp_to_total_num_patches)
    for mvpp, tot_patch_count in zip((1000, 10000), (273, 68)):
        edge_weights_all = extract(subject_id_list, example_dir,
                                   base_feature=base_feature,
                                   weight_method_list=[wm,],
                                   atlas='fsaverage', node_size=mvpp,
                                   smoothing_param=fwhm, out_dir=out_dir,
                                   return_results=True, num_procs=1)

        num_combinations = len(list(edge_weights_all))

        if num_combinations != len(subject_id_list):
            raise ValueError('mvpp: invalid count : # subjects')

        num_links = tot_patch_count * (tot_patch_count - 1) / 2
        for sub in subject_id_list:
            if edge_weights_all[(wm, sub)].size != num_links:
                raise ValueError('mvpp: invalid count : # links')


# test_multi_edge()
# test_multi_edge_CLI()
# test_empty_subject_list()
# test_run_no_IO()
# test_run_roi_stats_via_API()
# test_run_roi_stats_via_CLI()
# test_CLI_only_weight_or_stats()
# test_run_API_on_original_features()

test_atlas_parcel_subdivision()