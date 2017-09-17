import os
import sys
import stat
from os.path import join as pjoin, exists as pexists, realpath
import numpy as np
import pickle

from neuropredict import run_workflow
from pyradigm import MLDataset

#---------------------------------
# CHANGE THESE FOR YOUR PROCESSING
#---------------------------------

base_dir = '/u1/work/hpc3194'
dataset_list = ['PPMI', '4RTNI']

# freesurfer_dir = pjoin(proc_dir, 'freesurfer')
# subject_id_list = pjoin(target_list_dir, 'graynet.compute.list')

# list_of_datasets = [ '4RTNI', 'PPMI', 'ADNI' ]
# list_of_subject_lists = ['graynet.compute.list']*3

numeric_labels = {'CN' : 1, 'PARK' : 2, 'CBS' : 3, 'PSP': 4}

base_feature_list = [ 'freesurfer_thickness', 'freesurfer_curv' ] #  'freesurfer_thickness'
atlas = 'GLASSER2016' # 'FSAVERAGE' # 'GLASSER2016' #
fwhm = 10
node_size = None
feature_dim_expected = 360 # *(360-1)/2

roi_statistic = 'median'
compute_roi_statistic = True
expt_name = 'roistats'

num_splits_samples = 24 # 10.0
num_splits_stats = 1

roi_stat_list = np.array(['median', 'mode', 'mean', 'std', 'gmean', 'hmean',
                           'variation', 'entropy', 'skew', 'kurtosis'])

incomplete_processing, comb_nan_values = dict(), dict()

for base_feature in base_feature_list:
    print(' Processing {}'.format(base_feature))
    incomplete_processing[base_feature] = dict()
    comb_nan_values[base_feature] = dict()

    for stat_method in roi_stat_list:
        print('Gathering data for {}'.format(stat_method))
        expt_id = '{}_{}_{}_smoothing{}_size{}'.format(stat_method, base_feature, atlas, fwhm, node_size)

        flag_nan_exists = False
        flag_incomplete = False
        flag_unexpected = False
        dataset = MLDataset()

        incomplete_processing[base_feature][stat_method] = dict()
        comb_nan_values[base_feature][stat_method] = dict()
        for ds_name in dataset_list:
            print(' working on {}'.format(ds_name))
            proc_dir = pjoin(base_dir, ds_name, 'processed')
            out_dir = pjoin(proc_dir, 'graynet', '{}_{}_fwhm{}'.format(base_feature, atlas, fwhm))

            meta_list = pjoin(proc_dir, 'target_lists', 'meta_{}.csv'.format(ds_name))
            sample_ids, classes = run_workflow.get_metadata(meta_list)

            incomplete_processing[base_feature][stat_method][ds_name] = list()
            comb_nan_values[base_feature][stat_method][ds_name] = list()
            for sample in sample_ids:
                # 4139_bl_PPMI/median_freesurfer_thickness_GLASSER2016_smoothing10_sizeNone_roi_stats.csv
                feat_path = pjoin(out_dir, sample, '{}_roi_stats.csv'.format(expt_id))
                if pexists(feat_path):
                    data = np.genfromtxt(feat_path)
                    idx_nan = np.logical_not(np.isfinite(data))
                    local_flag_nan_exists = np.count_nonzero(idx_nan) > 0
                    if local_flag_nan_exists:
                        flag_nan_exists = True
                        comb_nan_values[base_feature][stat_method][ds_name].append(sample)
                        # print('NaNs found for {} {} {}'.format(ds_name, stat_method, sample))
                    elif len(data) >= feature_dim_expected:
                        dataset.add_sample(sample, data, numeric_labels[classes[sample]], class_id=classes[sample])
                    else:
                        flag_unexpected = True
                        incomplete_processing[base_feature][stat_method][ds_name].append(sample)
                else:
                    flag_incomplete = True
                    incomplete_processing[base_feature][stat_method][ds_name].append(sample)
                    # print('processing incomplete for {} {} {}'.format(ds_name, stat_method, sample))

        if flag_nan_exists or flag_incomplete or flag_unexpected:
            pass
            # print('{:20} {:25} - processing unusable; totally skipping it.'.format(base_feature, stat_method))
        else:
            print('{:20} {} \t- fully usable.'.format(base_feature, stat_method))
            dataset.description = '{}'.format(stat_method)
            out_path = pjoin(out_dir,'{}_roi_stats_{}.MLDataset.pkl'.format(base_feature,stat_method))
            dataset.save(out_path)

    # saving
    with open(pjoin(out_dir, 'incomplete_unusable_processing.pkl'), 'wb') as ipf:
        pickle.dump([incomplete_processing, comb_nan_values], ipf)

# reading
with open(pjoin(out_dir, 'incomplete_unusable_processing.pkl'), 'rb') as ipf:
    incomplete_processing, comb_nan_values = pickle.load(ipf)

# results
for base_feature in base_feature_list:
    for ds_name in dataset_list:
        for stat_method in roi_stat_list:
            print('{:20} {:5} {:25} {:5} {:5}'.format(base_feature, ds_name, stat_method,
                                                      len(incomplete_processing[base_feature][stat_method][ds_name]),
                                                      len(comb_nan_values[base_feature][stat_method][ds_name])))
