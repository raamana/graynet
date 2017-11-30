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

"""  # ---------------------------------#---------------------------------#---------------------------------

# ---------------------------------
# CHANGE THESE FOR YOUR PROCESSING
# ---------------------------------

dataset_name = '4RTNI' # 'ABIDE' # 'PPMI' # 'ADNI'  #   A short string identifying the larger dataset at a play

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
num_splits_samples = 10 # 25  # 10.0
num_splits_weights = 2

queue = 'abaqus.q'
mem = '8G'
walltime_per_subject = 0.5  # in hours -> 15 mins
cli_name = 'graynet'
job_type = 'multiedge_{}'.format(cli_name)
# ---------------------------------#---------------------------------#---------------------------------

atlas = 'fsaverage'  # 'glasser2016' # Choose one of 'fsaverage' and 'glasser2016'
fwhm = 5 # 10
num_bins = 25

base_feature_list = ['freesurfer_thickness', 'freesurfer_curv', 'freesurfer_sulc', 'freesurfer_area']
edge_range_dict = {'freesurfer_thickness': (0.0, 5.0),
                   'freesurfer_curv'     : (-0.3, +0.3),
                   'freesurfer_sulc'     : (-1.5, +1.5),
                   'freesurfer_area'     : (0.0, 1.5)
                   }

feat_names_flat = ' '.join(base_feature_list)
multi_edge_range_flat = ''
for feat in base_feature_list:
    multi_edge_range_flat += ' {} {}'.format(edge_range_dict[feat][0], edge_range_dict[feat][1])

expt_prefix = 'thk_curv_sulc_area'
# You can choose only one or multiple, but keep them enclosed as a list or array.
weights = np.array(['chebyshev', 'chi_square', 'correlate', 'cosine', 'euclidean',
                           'histogram_intersection', 'jensen_shannon', 'manhattan', 'minowski', 'relative_deviation'])

summary_stat_list = [ 'prod', 'median', 'amax', 'amin', 'gmean', 'std' ]
summary_stat = ' '.join(summary_stat_list)

cluster_type = 'SGE'
pe_name = 'dist.pe'  # distributed parallel env, not shared mem


def specify_hpc_resources(mem, queue, num_procs, job_dir, job_log, cluster_type='SGE'):
    "returns lines to include in job scripts to specify resources"

    lines = list()
    if cluster_type.upper() in ['SGE', 'SUNGRIDENGINE']:
        parallel_str = ''
        if num_procs > 1:
            parallel_str = '-pe {} {}'.format(pe_name, num_procs)
        lines.append('#$ -l mf={} -q {} {}'.format(mem, queue, parallel_str))
        lines.append('#$ -wd {}'.format(job_dir))
        lines.append('#$ -j yes -o {}'.format(job_log))
        lines.append('cd {}\n'.format(job_dir))
    else:
        # add other lines below accommodating to your own HPC
        # contact me if you need help - its not an obscure HPC, i can quickly help you with this.
        raise TypeError

    return '\n'.join(lines)


# ---------------------------------
# END MAKING CHANGES
# ---------------------------------


out_dir = pjoin(proc_dir, 'graynet', '{}_{}_fwhm{}_nbins{}'.format(expt_prefix, atlas, fwhm, num_bins))
os.makedirs(out_dir, exist_ok=True)
print('Saving the jobs to\n{}'.format(out_dir))


def make_dirs(dir_list):
    "Make multiple directories at once"

    for dir in dir_list:
        if not pexists(dir):
            os.makedirs(dir)


def make_cli_call(cli_name, subject_id_list, base_feature, freesurfer_dir,
                  weight_method, num_bins, edge_range_flat, summary_stat, atlas, fwhm, out_proc_dir, num_procs=2):
    cmd = '{} -s {} -f {} -i {} -w {} ' \
          '--do_multi_edge --multi_edge_range {} --summary_stat {} ' \
          '-b {} -a {} -p {} -o {} -c {} \n'.format(cli_name,
                                                    realpath(subject_id_list),
                                                    base_feature,
                                                    realpath(freesurfer_dir),
                                                    weight_method, edge_range_flat, summary_stat,
                                                    num_bins, atlas,
                                                    fwhm, out_proc_dir, num_procs)

    return cmd


def make_job(subject_id_list, freesurfer_dir,
             base_feature, weight_method, num_bins, edge_range, summary_stat,
             atlas, fwhm, out_proc_dir, job_dir, job_name, num_procs):
    "Creates graynet job for running on HPC"

    str_list_weight_method = ' '.join(weight_method)

    job_file = pjoin(job_dir, '{}.{}.job'.format(job_name, job_type))
    job_log = pjoin(job_dir, '{}.{}.log'.format(job_name, job_type))
    if pexists(job_file):
        os.remove(job_file)
    with open(job_file, 'w') as jf:
        jf.write('#!/bin/bash\n')
        jf.write(specify_hpc_resources(mem, queue, num_procs, job_dir, job_log))
        jf.write(make_cli_call(cli_name, realpath(subject_id_list), base_feature, realpath(freesurfer_dir),
                               str_list_weight_method, num_bins, edge_range, summary_stat, atlas, fwhm, realpath(out_proc_dir),
                               num_procs))

    st = os.stat(job_file)
    os.chmod(job_file, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    return job_file


job_dir = pjoin(out_dir, 'PBS')
make_dirs([job_dir, ])

# astype(str) is important to avoid byte strings
id_list = np.atleast_1d(np.genfromtxt(subject_id_list, dtype=str).astype(str))
num_samples = len(id_list)
num_weights = len(weights)

print('{} samples in {} splits\n{} weights in {} splits'.format(num_samples, num_splits_samples, num_weights, num_splits_weights))

subsets_sub_id = np.array_split(id_list,num_splits_samples)
subsets_weight = np.array_split(weights,num_splits_weights)

for ww in range(num_splits_weights):

    for ss in range(num_splits_samples):

        subset_list_path = pjoin(job_dir, '{}_split{}_samples.txt'.format(dataset_name.lower(), ss))
        with open(subset_list_path, 'w') as sf:
            sf.write('\n'.join(subsets_sub_id[ss]))

        job_name = 'split_w{}_s{}'.format(ww, ss)
        job_file = make_job(subset_list_path, freesurfer_dir,
                            feat_names_flat, subsets_weight[ww], num_bins,
                            multi_edge_range_flat, summary_stat,
                            atlas, fwhm, out_dir, job_dir, job_name, num_procs)

        print('generated jobs for {}'.format(job_name))

print('Job scripts have been generated in \n {}'.format(job_dir))
print('\nPlease change to that directory and submit them using ./qsubMany *.graynet.job')
print('If you are not on Linux or not using SGE/Torque, you may have to adapt the qsubMany script also.')
