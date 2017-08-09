
import os
from os.path import join as pjoin
import sys
import socket
import numpy as np
import nibabel as nib


if socket.gethostname().startswith('SQuark'):
    fs_subjects_dir = '/mnt/opt/freesurfer/subjects'
else:
    fs_subjects_dir = '/mnt/opt/freesurfer/subjects'

atlas_name = 'fsaverage'
atlas_dir = pjoin(fs_subjects_dir, atlas_name)

hemi_list = ['lh', 'rh']
atlas, coords, faces, info = dict(), dict(), dict(), dict()

for hemi in hemi_list:
    hemi_path = os.path.join(atlas_dir,'surf','lh.orig')
    coords[hemi], faces[hemi], info[hemi] = nib.freesurfer.io.read_geometry(hemi_path, read_metadata=True)

num_vertices_left_hemi = coords['lh'].shape[0]

coords['whole'] = np.vstack((coords['lh'], coords['rh']))
faces['whole'] = np.vstack((faces['lh'], faces['rh'] + num_vertices_left_hemi))

# labels, ctab, names = nib.freesurfer.io.read_annot(pjoin(mmp_fs,'lh.HCP-MMP1.annot'))

# reading annotations
annot = dict()
for hemi in hemi_list:
    annot[hemi] = dict()
    annot_path = pjoin(atlas_dir, 'label', '{}.aparc.annot'.format(hemi))
    annot[hemi]['labels'], annot[hemi]['ctab'], annot[hemi]['names'] = nib.freesurfer.io.read_annot(annot_path, orig_ids=True)

labels_to_remove = ['corpuscallosum', 'unknown']
null_label= 0

def ismember(A, B):
    B_unique_sorted, B_idx = np.unique(B, return_index=True)
    B_in_A_bool = np.in1d(B_unique_sorted, A, assume_unique=True)

for hemi in hemi_list:
    cortex_label_path = pjoin(atlas_dir, 'label', '{}.cortex.label'.format(hemi))
    cortex_label = nib.freesurfer.io.read_label(cortex_label_path)

    mask_for_cortex = np.in1d(annot[hemi]['labels'], cortex_label, assume_unique=True)


print('lets do something')