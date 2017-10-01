#!/usr/bin/env python

import os
import sys
import stat
from os.path import join as pjoin, exists as pexists, realpath
import numpy as np

"""#---------------------------------#---------------------------------#---------------------------------

The following script helps you generating job files for processing the large computation on a cluster.

Just input details regarding 
- where your Freesurfer processing folder (SUBJECTS_DIR) is 
- and what subjects to process
- what feature is your "base features" to extract pair-wise network level features from

We expect Freesurfer processing to have been done with -qcache flag with an FWHM value of 10. 
The Qcache recon-all flag completes the following:
    - resample data (e.g. thickness, curv, sulc) onto the average subject (called fsaverage)
    - smooth it at a range of FWHM (full-width/half-max) values, usually 0, 5, 10, 15, 20, and 25mm. 
Check https://surfer.nmr.mgh.harvard.edu/fswiki/qcache for more info

graynet comes packaged with parcellations (info on what vertices belong to what ROIs) of fsaverage and Glasser2016. So you simply need to select which ones you would like. 

There are plans to 
 - integrate arbitrary atlases - for example, dataset- and age-specific atlases, as long as their info is organized in a Freesurfer structure.
 - offer subdivision the defined ROIs to control their size.
However, this is still a work in progress. Please email me and I will be happy to make it happen sooner than later. 

"""#---------------------------------#---------------------------------#---------------------------------

#---------------------------------
# CHANGE THESE FOR YOUR PROCESSING
#---------------------------------

dataset_name = 'ADNI' # 'PPMI' #  '4RTNI' # A short string identifying the larger dataset at a play

# the following paths can be derived from any number of ways -
#   just ensure freesurfer_dir and subject_id_list are defined and exist.
base_dir = '/u1/work/hpc3194'
proc_dir = pjoin(base_dir, dataset_name, 'processed')
target_list_dir = pjoin(proc_dir, 'target_lists')

# path to SUBJECTS_DIR
freesurfer_dir = pjoin(proc_dir, 'freesurfer')
# this must be file with one subject ID per line
subject_id_list = pjoin(target_list_dir, 'graynet.compute.list')

num_procs = 4
num_splits_samples = 25.0 # 10.0
num_splits_weights =  2.0

queue = 'abaqus.q'
mem='8G'
walltime_per_subject = 0.5 # in hours -> 15 mins
cli_name = 'graynet'

#---------------------------------#---------------------------------#---------------------------------

base_feature = 'freesurfer_thickness' # Choose one of 'freesurfer_thickness', 'freesurfer_curv', 'freesurfer_sulc'
atlas = 'FSAVERAGE' #  'GLASSER2016' # Choose one of 'FSAVERAGE' and 'GLASSER2016'
fwhm = 10

# number of bins used in constructing a histograms used to compute pair-wise distances
num_bins = 25

# this is to control the full range of values within which histogram bins are defined
# I think these ranges are reasonable for thickness and curvature based on what I have seen in my datasets,
# but feel free to change using your own measurements and datasets.
edge_range_predefined = {'freesurfer_thickness': (0, 5), 'freesurfer_curv': (-0.3, +0.3)}
edge_range = edge_range_predefined[base_feature]

# You can choose only one or multiple, but keep them enclosed as a list or array.
histogram_dist = np.array(['chebyshev', 'chi_square', 'correlate', 'cosine', 'euclidean',
                           'histogram_intersection', 'jensen_shannon', 'manhattan', 'minowski',  'relative_deviation'])


cluster_type = 'SGE'
pe_name = 'dist.pe' # distributed parallel env, not shared mem

def specify_hpc_resources(mem, queue, num_procs, job_dir, job_log, cluster_type='SGE'):
    "returns lines to include in job scripts to specify resources"

    lines = list()
    if cluster_type.upper() in ['SGE', 'SUNGRIDENGINE']:
        parallel_str=''
        if num_procs > 1:
                parallel_str='-pe {} {}'.format(pe_name,num_procs)
        lines.append('#$ -l mf={} -q {} {} -wd {} -j yes -o {}'.format(mem, queue, parallel_str, job_dir, job_log))
        lines.append('cd {}\n'.format(job_dir))
    else:
        # add other lines below accommodating to your own HPC
        # contact me if you need help - its not an obscure HPC, i can quickly help you with this.
        raise TypeError

    return '\n'.join(lines)


#---------------------------------
# END MAKING CHANGES
#---------------------------------


out_dir = pjoin(proc_dir, 'graynet', '{}_{}_fwhm{}_range{}_{}_nbins{}'.format(base_feature, atlas, fwhm, edge_range[0], edge_range[1], num_bins))


if not pexists(out_dir):
    os.makedirs(out_dir)
print('Saving the jobs to\n{}'.format(out_dir))

def make_dirs(dir_list):
    "Make multiple directories at once"

    for dir in dir_list:
        if not pexists(dir):
            os.makedirs(dir)

def make_cli_call(cli_name, subject_id_list, base_feature, freesurfer_dir, 
                weight_method, num_bins, edge_range, atlas, fwhm, out_proc_dir, num_procs=2):

    cmd = '{} -s {} -f {} -i {} -w {} -e {} {} -b {} -a {} -p {} -o {} -c {} \n'.format(cli_name,
            realpath(subject_id_list), base_feature, realpath(freesurfer_dir),
            weight_method, edge_range[0], edge_range[1], num_bins, atlas, fwhm, out_proc_dir, num_procs)

    return cmd


def make_job(subject_id_list, freesurfer_dir,
             base_feature, weight_method, num_bins, edge_range,
             atlas, fwhm, out_proc_dir, job_dir, job_name, num_procs):
    "Creates graynet job for running on HPC"

    str_list_weight_method = ' '.join(weight_method)

    job_file = pjoin(job_dir, '{}.graynet.job'.format(job_name))
    job_log  = pjoin(job_dir, '{}.graynet.log'.format(job_name))
    if pexists(job_file):
        os.remove(job_file)
    with open(job_file, 'w') as jf:
        jf.write('#!/bin/bash\n')
        jf.write(specify_hpc_resources(mem, queue, num_procs, job_dir, job_log))
        jf.write(make_cli_call(cli_name,realpath(subject_id_list), base_feature, realpath(freesurfer_dir),
            str_list_weight_method, num_bins, edge_range, atlas, fwhm, realpath(out_proc_dir), num_procs))

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
                            atlas, fwhm, out_dir, job_dir, job_name, num_procs)

        print('generated jobs for {}'.format(job_name))

print('Job scripts have been generated in \n {}'.format(job_dir))
print('\nPlease change to that directory and submit them using ./qsubMany *.graynet.job')
print('If you are not on Linux or not using SGE/Torque, you may have to adapt the qsubMany script also.')
