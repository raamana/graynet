
import os
import sys
import stat
from os.path import join as pjoin, exists as pexists, realpath
import numpy as np
from pyradigm import MLDataset
from graynet.run_workflow import import_features
from neuropredict.run_workflow import get_metadata

base_dir = '/u1/work/hpc3194'
dataset_name = '4RTNI' # 'PPMI' #

# list_of_datasets = [ '4RTNI', 'PPMI', 'ADNI' ]
# list_of_subject_lists = ['graynet.compute.list']*3

proc_dir = pjoin(base_dir, dataset_name, 'processed')
freesurfer_dir = pjoin(proc_dir, 'freesurfer')
target_list_dir = pjoin(proc_dir, 'target_lists')

vis_out_dir = pjoin(proc_dir, 'visualizations')
if not pexists(vis_out_dir):
    os.mkdir(vis_out_dir)

subject_id_list = pjoin(target_list_dir, 'graynet.compute.list')
# meta_file = pjoin(target_list_dir, 'meta_4RTNI_and_PPMI.csv')
meta_file = pjoin(target_list_dir, 'meta_{}.csv'.format(dataset_name))

base_feature = 'freesurfer_thickness' # 'freesurfer_curv' #  'freesurfer_thickness'
atlas = 'GLASSER2016' # 'FSAVERAGE' # 'GLASSER2016' #
fwhm = 10

# id_list = np.atleast_1d(np.genfromtxt(subject_id_list, dtype=str).astype(str))
# num_samples = len(id_list)

id_list, classes = get_metadata(meta_file)
class_set = list(set(classes.values()))
class_set.sort()
labels = { sub : class_set.index(cls) for sub, cls in classes.items() }

id_data = import_features(freesurfer_dir, id_list, base_feature)

ds = MLDataset(data=id_data, labels=labels, classes=classes)

out_path = pjoin(vis_out_dir, 'raw_features_{}_{}.MLDataset.pkl'.format(base_feature, '_'.join(class_set)))
ds.save(out_path)

data, lbl, ids = ds.data_and_labels()
print('{} {}\n   min: {:.3f}\n   max: {:.3f}\n 0.1% : {:.3f}\n99.9% : {:.3f}'.format(dataset_name, base_feature, np.min(data), np.max(data),
                                                    np.percentile(data, 0.1), np.percentile(data, 99.9)))