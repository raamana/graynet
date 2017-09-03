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
dataset_name = '4RTNI' #

# list_of_datasets = [ '4RTNI', 'PPMI', 'ADNI' ]
# list_of_subject_lists = ['graynet.compute.list']*3

proc_dir = pjoin(base_dir, dataset_name, 'processed')
freesurfer_dir = pjoin(proc_dir, 'freesurfer')
target_list_dir = pjoin(proc_dir, 'target_lists')

subject_id_list = pjoin(target_list_dir, 'graynet.compute.list')
print('test run with 5 subjects only')

base_feature = 'freesurfer_thickness'
atlas = 'GLASSER2016' # 'FSAVERAGE' # 'GLASSER2016' #
fwhm = 10

out_dir = pjoin(proc_dir, 'graynet', '{}_{}_fwhm{}'.format(base_feature, atlas, fwhm))

cluster_type = 'SGE'

def specify_hpc_resources(mem, queue, job_dir, cluster_type='SGE'):
    "returns lines to include in job scripts to specify resources"

    lines = list()
    if cluster_type.upper() in ['SGE', 'SUNGRIDENGINE']:
        lines.append('#$ -l mf={} -q {}'.format(mem, queue))
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

def make_dirs(dir_list):
    "Make multiple directories at once"

    for dir in dir_list:
        if not pexists(dir):
            os.makedirs(dir)

def make_cli_call(cli_name, subject_id_list, base_feature, freesurfer_dir,
            weight_method, atlas, fwhm, out_proc_dir):

    cmd = '{} -s {} -f {} -i {} -w {} -a {} -p {} -o {}\n'.format(cli_name,
            realpath(subject_id_list), base_feature, realpath(freesurfer_dir),
            weight_method, atlas, fwhm, out_proc_dir)

    return cmd


def make_job(subject_id_list, freesurfer_dir,
             base_feature, weight_method,
             atlas, fwhm,
             out_proc_dir, job_dir, job_name):
    "Creates graynet job for running on HPC"

    queue = 'abaqus.q'
    mem='4G'
    walltime_per_subject = 0.25 # in hours -> 15 mins
    cli_name = 'graynet'

    sub_list = np.atleast_1d(np.genfromtxt(subject_id_list))
    num_subjects = len(sub_list)
    walltime = int(np.round(walltime_per_subject*num_subjects, 0))

    str_list_weight_method = ' '.join(weight_method)

    job_file = pjoin(job_dir, '{}.job'.format(job_name))
    if pexists(job_file):
        os.remove(job_file)
    with open(job_file, 'w') as jf:
        jf.write('#!/bin/bash\n')
        jf.write(specify_hpc_resources(mem, queue, job_dir, cluster_type))
        jf.write(make_cli_call(cli_name,realpath(subject_id_list), base_feature, realpath(freesurfer_dir),
            str_list_weight_method, atlas, fwhm, realpath(out_proc_dir)))

    st = os.stat(job_file)
    os.chmod(job_file, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    return job_file


histogram_dist = np.array([
    'chebyshev', 'chebyshev_neg', 'chi_square',
    'correlate', 'correlate_1',
    'cosine', 'cosine_1', 'cosine_2', 'cosine_alt',
    'euclidean', 'fidelity_based',
    'histogram_intersection', 'histogram_intersection_1',
    'jensen_shannon', 'kullback_leibler', 'manhattan', 'minowski',
    'noelle_1', 'noelle_2', 'noelle_3', 'noelle_4', 'noelle_5',
    'relative_bin_deviation', 'relative_deviation'])

num_splits_samples = 3.0 # 10.0
num_splits_weights = 3.0

job_dir = pjoin(out_dir, 'PBS')
make_dirs([job_dir, ])

num_weights = len(histogram_dist)
id_list = np.atleast_1d(np.loadtxt(subject_id_list, dtype=str))
num_samples = len(id_list)

num_samples_per_job = max(1,np.int64(np.ceil(num_samples/num_splits_samples)))
num_weights_per_job = max(1,np.int64(np.ceil(num_weights/num_splits_weights)))

wt_count = 0
for ww in range(int(num_splits_weights)):
    subset_weights = histogram_dist[wt_count:min(wt_count+num_weights_per_job, num_weights)]
    wt_count = wt_count + num_weights_per_job

    print(' {} \n {}'.format(ww, subset_weights))

    sub_count = 0
    for ss in range(int(num_splits_samples)):
        subset_samples = id_list[sub_count:min(sub_count+num_samples_per_job,num_samples)]
        sub_count = sub_count + num_samples_per_job

        print(' {} \n {}'.format(ss, subset_samples))

        subset_list_path = pjoin(out_dir,'split{}_samples.txt'.format(ss))
        with open(subset_list_path, 'w') as sf:
            sf.write('\n'.join(subset_samples))

        job_name = 'split{}{}_{}'.format(ww,ss,'_'.join(subset_weights))
        job_file = make_job(subset_list_path, freesurfer_dir,
                            base_feature, subset_weights,
                            atlas, fwhm, out_dir, job_dir, job_name)

        print('generated jobs for {}'.format(job_name))