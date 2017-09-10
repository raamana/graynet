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
dataset_name = '4RTNI' # 'PPMI' #

# freesurfer_dir = pjoin(proc_dir, 'freesurfer')


# subject_id_list = pjoin(target_list_dir, 'graynet.compute.list')

base_feature = 'freesurfer_thickness' # 'freesurfer_curv' #  'freesurfer_thickness'
atlas = 'GLASSER2016' # 'FSAVERAGE' # 'GLASSER2016' #
fwhm = 10
node_size = None

histogram_dist = np.array([
    'chebyshev', 'chebyshev_neg', 'chi_square',
    'correlate', 'correlate_1',
    'cosine', 'cosine_1', 'cosine_2', 'cosine_alt',
    'euclidean', 'fidelity_based',
    'histogram_intersection', 'histogram_intersection_1',
    'jensen_shannon', 'kullback_leibler', 'manhattan', 'minowski',
    'noelle_1', 'noelle_2', 'noelle_3', 'noelle_4', 'noelle_5',
    'relative_bin_deviation', 'relative_deviation'])

incomplete_processing, comb_nan_values = dict(), dict()

for weight_method in histogram_dist: # ['chi_square', 'cosine', 'euclidean', 'histogram_intersection']:
    # print('Gathering data for {}'.format(weight_method))

    expt_id = '{}_{}_smoothing{}_size{}_edgeweight_{}'.format(base_feature, atlas, fwhm, node_size,
                                                              weight_method)

    flag_nan_exists = False
    flag_incomplete = False
    dataset = MLDataset()

    incomplete_processing[weight_method] = dict()
    comb_nan_values[weight_method] = dict()
    for ds_name in dataset_list:
        proc_dir = pjoin(base_dir, dataset_name, 'processed')
        out_dir = pjoin(proc_dir, 'graynet', '{}_{}_fwhm{}'.format(base_feature, atlas, fwhm))

        meta_list = pjoin(proc_dir, 'target_lists', 'meta_{}.csv'.format(ds_name))
        sample_ids, classes = run_workflow.get_metadata(meta_list)

        incomplete_processing[weight_method][ds_name] = list()
        comb_nan_values[weight_method][ds_name] = list()
        for sample in sample_ids:
            feat_path = pjoin(out_dir, sample, '{}_graynet.csv'.format(expt_id))
            if pexists(feat_path):
                data = np.genfromtxt(feat_path)
                idx_nan = np.isnan(data)
                local_flag_nan_exists = np.count_nonzero(idx_nan) > 0
                if local_flag_nan_exists:
                    flag_nan_exists = True
                    comb_nan_values[weight_method][ds_name].append(sample)
                    # print('NaNs found for {} {} {}'.format(ds_name, weight_method, sample))
                else:
                    dataset.add_sample(sample, data, classes[sample], class_id=classes[sample])
            else:
                flag_incomplete = True
                # print('processing incomplete for {} {} {}'.format(ds_name, weight_method, sample))
                incomplete_processing[weight_method][ds_name].append(sample)

    if flag_nan_exists or flag_incomplete:
        print('NaNs/processing incmplete for {} - totally skipping it.'.format(weight_method))
    else:
        print(' {} {} - fully usable.'.format(ds_name, weight_method))
        dataset.description = 'weight {}'.format(weight_method)
        out_path = pjoin(out_dir,'{}.MLDataset.pkl'.format(weight_method))
        dataset.save(out_path)

# saving
with open(pjoin('.', 'incomplete_processing.pkl'), 'wb') as ipf:
    pickle.dump([incomplete_processing, comb_nan_values], ipf)

# reading
with open(pjoin('.', 'incomplete_processing.pkl'), 'rb') as ipf:
    incomplete_processing, comb_nan_values = pickle.load(ipf)

for ds_name in dataset_list:
    for wm in histogram_dist:
        print('{} {:25} {:5} {:5}'.format(ds_name, wm,
            len(incomplete[weight_method][ds_name]), len(comb_nan[weight_method][ds_name])))
