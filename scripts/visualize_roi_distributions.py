
import os
import sys
import stat
from os.path import join as pjoin, exists as pexists, realpath
import numpy as np
from pyradigm import MLDataset
from graynet.run_workflow import import_features
from neuropredict.run_workflow import get_metadata
import traceback

base_dir = '/u1/work/hpc3194'

features_freesurfer = ['freesurfer_thickness',
                       'freesurfer_sulc',
                       'freesurfer_curv',
                       'freesurfer_area',]

atlas = 'fsaverage' # 'GLASSER2016' # 'FSAVERAGE' # 'GLASSER2016' #
fwhm = 10

for dataset_name in ['4RTNI', 'PPMI' ]:

    proc_dir = pjoin(base_dir, dataset_name, 'processed')
    freesurfer_dir = pjoin(proc_dir, 'freesurfer')
    target_list_dir = pjoin(proc_dir, 'target_lists')

    vis_out_dir = pjoin(proc_dir, 'visualizations')
    if not pexists(vis_out_dir):
        os.mkdir(vis_out_dir)

    subject_id_list = pjoin(target_list_dir, 'graynet.compute.list')
    # meta_file = pjoin(target_list_dir, 'meta_4RTNI_and_PPMI.csv')
    meta_file = pjoin(target_list_dir, 'meta_{}.csv'.format(dataset_name))

    for base_feature in features_freesurfer:

        id_list, classes = get_metadata(meta_file)
        class_set = list(set(classes.values()))
        class_set.sort()
        labels = { sub : class_set.index(cls) for sub, cls in classes.items() }

        out_path = pjoin(vis_out_dir, 'raw_features_{}_{}.MLDataset.pkl'.format(base_feature, '_'.join(class_set)))

        try:
            ds = MLDataset(filepath=out_path)
        except:
            traceback.print_exc()
            id_data = import_features(freesurfer_dir, id_list, base_feature)
            ds = MLDataset(data=id_data, labels=labels, classes=classes)
            ds.save(out_path)


        data, lbl, ids = ds.data_and_labels()
        print('{} {}\n min : {:.4f}\n max : {:.4f}'.format(dataset_name, base_feature, np.min(data), np.max(data)))
        for perc in [ 1, 5, 95, 99 ] :
            print('{:3d}% : {:10.4f}'.format(perc, np.percentile(data, perc)))
