
import os
from os.path import join as pjoin, exists as pexists
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

def get_atlas_annot(atlas_name=None):
    "High level wrapper to get all the info just by using a name."

    if atlas_name in [None, 'None']:
        atlas_name = 'GLASSER2016'

    atlas_name = atlas_name.upper()

    if atlas_name in ['GLASSER2016']:
        atlas_path = '/data1/strother_lab/praamana/glasser2016/fsaverage_annot_figshare3498446'
    elif atlas_name in ['FSAVERAGE']:
        atlas_path = '/mnt/opt/freesurfer/subjects'
        if not pexists(atlas_path):
            atlas_path = '/opt/freesurfer/subjects'
    else:
        raise NotImplementedError('Requested atlas is not implemented or unreadable.')

    annot = read_atlas_annot(atlas_path)
    return annot


def read_atlas_annot(atlas_dir, hemi_list=None):
    " Returns atlas annotations "

    if hemi_list is None:
        hemi_list = ['lh', 'rh']

    annot = dict()
    for hemi in hemi_list:
        annot[hemi] = dict()
        annot_path = pjoin(atlas_dir, 'label', '{}.aparc.annot'.format(hemi))
        annot[hemi]['labels'], annot[hemi]['ctab'], annot[hemi]['names'] = nib.freesurfer.io.read_annot(annot_path, orig_ids=True)

    return annot


def read_atlas(atlas_dir, hemi_list=None):
    " Script to read the pre-computed parcellations for fsaverage and HCP-MMP-1.0 "

    if hemi_list is None:
        hemi_list = ['lh', 'rh']

    coords, faces, info = dict(), dict(), dict()

    for hemi in hemi_list:
        hemi_path = os.path.join(atlas_dir,'surf','lh.orig')
        coords[hemi], faces[hemi], info[hemi] = nib.freesurfer.io.read_geometry(hemi_path, read_metadata=True)

    num_vertices_left_hemi = coords['lh'].shape[0]

    coords['whole'] = np.vstack((coords['lh'], coords['rh']))
    faces['whole'] = np.vstack((faces['lh'], faces['rh'] + num_vertices_left_hemi))

    # labels, ctab, names = nib.freesurfer.io.read_annot(pjoin(mmp_fs,'lh.HCP-MMP1.annot'))

    annot = read_atlas_annot(atlas_dir, hemi_list)

    return coords, faces, annot


def subdivide_cortex(atlas_dir):
    "Subdivides the given cortical parcellation (each label into smaller patches)"

    raise NotImplementedError('This function has not been implemented yet.')

    # noinspection PyUnreachableCode
    coords, faces, annot = read_atlas(atlas_dir)

    labels_to_remove = ['corpuscallosum', 'unknown']
    null_label= 0

    def ismember(A, B):
        B_unique_sorted, B_idx = np.unique(B, return_index=True)
        B_in_A_bool = np.in1d(B_unique_sorted, A, assume_unique=True)

    cortex_label = dict()
    for hemi in hemi_list:
        cortex_label_path = pjoin(atlas_dir, 'label', '{}.cortex.label'.format(hemi))
        cortex_label[hemi] = nib.freesurfer.io.read_label(cortex_label_path)

    # # cortex_label[hemi] is an index into annot[hemi]['labels']

        mask_for_cortex = np.in1d(annot[hemi]['labels'], cortex_label, assume_unique=True)


print('lets do something')