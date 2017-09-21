#!/usr/bin/env python

import os
import sys
import stat
from os.path import join as pjoin, exists as pexists, realpath
import numpy as np

#---------------------------------
# CHANGE THESE FOR YOUR PROCESSING
#---------------------------------

base_dir = '/u1/work/hpc3194'
dataset_name = '4RTNI' # 'PPMI' #

# list_of_datasets = [ '4RTNI', 'PPMI', 'ADNI' ]
# list_of_subject_lists = ['graynet.compute.list']*3

proc_dir = pjoin(base_dir, dataset_name, 'processed')
freesurfer_dir = pjoin(proc_dir, 'freesurfer')
target_list_dir = pjoin(proc_dir, 'target_lists')

subject_id_list = pjoin(target_list_dir, 'graynet.compute.list')

base_feature = 'freesurfer_thickness' # 'freesurfer_curv' #
atlas = 'GLASSER2016' # 'FSAVERAGE' # 'GLASSER2016' #
fwhm = 10

num_bins = 25
edge_range_predefined = {'freesurfer_thickness': (0, 5), 'freesurfer_curv': (-0.3, +0.3)}
edge_range = edge_range_predefined[base_feature]

out_dir = pjoin(proc_dir, 'graynet', '{}_{}_fwhm{}_range{}_{}_nbins{}'.format(base_feature, atlas, fwhm, edge_range[0], edge_range[1], num_bins))

histogram_dist = np.array(['chebyshev', 'chi_square', 'correlate', 'cosine', 'euclidean',
                           'histogram_intersection', 'jensen_shannon', 'manhattan', 'minowski',  'relative_deviation'])

# histogram_dist = np.array([
#     'chebyshev', 'chebyshev_neg', 'chi_square',
#     'correlate', 'correlate_1',
#     'cosine', 'cosine_1', 'cosine_2', 'cosine_alt',
#     'euclidean', 'fidelity_based',
#     'histogram_intersection', 'histogram_intersection_1',
#     'jensen_shannon', 'kullback_leibler', 'manhattan', 'minowski',
#     'noelle_1', 'noelle_2', 'noelle_3', 'noelle_4', 'noelle_5',
#     'relative_bin_deviation', 'relative_deviation'])

cluster_type = 'SGE'

def specify_hpc_resources(mem, queue, job_dir, job_log, cluster_type='SGE'):
    "returns lines to include in job scripts to specify resources"

    lines = list()
    if cluster_type.upper() in ['SGE', 'SUNGRIDENGINE']:
        lines.append('#$ -l mf={} -q {} -wd {} -j yes -o {}'.format(mem, queue, job_dir, job_log))
        lines.append('cd {}\n'.format(job_dir))
        # add other lines below accommodating to your own HPC
    else:
        raise TypeError

    return '\n'.join(lines)


#---------------------------------
# END MAKING CHANGES
#---------------------------------

if not pexists(out_dir):
    os.makedirs(out_dir)
print('Saving the jobs to\n{}'.format(out_dir))

def make_dirs(dir_list):
    "Make multiple directories at once"

    for dir in dir_list:
        if not pexists(dir):
            os.makedirs(dir)

def make_cli_call(cli_name, subject_id_list, base_feature, freesurfer_dir,
            weight_method, num_bins, edge_range, atlas, fwhm, out_proc_dir):

    cmd = '{} -s {} -f {} -i {} -w {} -e {} {} -b {} -a {} -p {} -o {}\n'.format(cli_name,
            realpath(subject_id_list), base_feature, realpath(freesurfer_dir),
            weight_method, edge_range[0], edge_range[1], num_bins, atlas, fwhm, out_proc_dir)

    return cmd


def make_job(subject_id_list, freesurfer_dir,
             base_feature, weight_method, num_bins, edge_range,
             atlas, fwhm,
             out_proc_dir, job_dir, job_name):
    "Creates graynet job for running on HPC"

    queue = 'abaqus.q'
    mem='4G'
    walltime_per_subject = 0.5 # in hours -> 15 mins
    cli_name = 'graynet'

    str_list_weight_method = ' '.join(weight_method)

    job_file = pjoin(job_dir, '{}.job'.format(job_name))
    job_log  = pjoin(job_dir, '{}.log'.format(job_name))
    if pexists(job_file):
        os.remove(job_file)
    with open(job_file, 'w') as jf:
        jf.write('#!/bin/bash\n')
        jf.write(specify_hpc_resources(mem, queue, job_dir, job_log))
        jf.write(make_cli_call(cli_name,realpath(subject_id_list), base_feature, realpath(freesurfer_dir),
            str_list_weight_method, num_bins, edge_range, atlas, fwhm, realpath(out_proc_dir)))

    st = os.stat(job_file)
    os.chmod(job_file, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    return job_file


job_dir = pjoin(out_dir, 'PBS')
make_dirs([job_dir, ])

num_weights = len(histogram_dist)
# astype(str) is important to avoid byte strings
id_list = np.atleast_1d(np.genfromtxt(subject_id_list, dtype=str).astype(str))
num_samples = len(id_list)

print('{} samples and {} weights'.format(num_samples, num_weights))

num_splits_samples = 50.0 # 10.0
num_splits_weights = 2.0

num_samples_per_job = max(1,np.int64(np.ceil(num_samples/num_splits_samples)))
num_weights_per_job = max(1,np.int64(np.ceil(num_weights/num_splits_weights)))

print('{} samples/job and {} weights/job'.format(num_samples_per_job, num_weights_per_job))

wt_count, end_idx1 = 0, 0
for ww in range(int(num_splits_weights)):
    if end_idx1 > num_weights:
        break
    end_idx1 = wt_count+num_weights_per_job
    subset_weights = histogram_dist[wt_count:end_idx1]
    wt_count = end_idx1

    sub_count, end_idx2 = 0, 0
    for ss in range(int(num_splits_samples)):
        if end_idx2 > num_samples:
            break
        end_idx2 = sub_count+num_samples_per_job+1
        subset_samples = id_list[sub_count:end_idx2]
        sub_count = end_idx2

        subset_list_path = pjoin(job_dir,'{}_{}th_split_samples.txt'.format(dataset_name.lower(),ss))
        with open(subset_list_path, 'w') as sf:
            sf.write('\n'.join(subset_samples))

        job_name = 'split_w{}_s{}_{}'.format(ww,ss,'_'.join(subset_weights))
        job_file = make_job(subset_list_path, freesurfer_dir,
                            base_feature, subset_weights, num_bins, edge_range,
                            atlas, fwhm, out_dir, job_dir, job_name)

        print('generated jobs for {}'.format(job_name))