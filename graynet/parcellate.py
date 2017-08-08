
import os
import sys
import numpy as np
import nibabel as nib

atlas_dir = '/mnt/opt/freesurfer/subjects/fsaverage'

left_hemi = os.path.join(atlas_dir,'surf','lh.orig')
left_coords, left_faces, left_info = nib.freesurfer.io.read_geometry(left_hemi, read_metadata=True)

right_hemi = os.path.join(atlas_dir,'surf','rh.orig')
right_coords, right_faces, right_info = nib.freesurfer.io.read_geometry(right_hemi, read_metadata=True)

num_vertices_left_hemi = left_coords.shape[0]
num_vertices_right_hemi = right_coords.shape[0]
assert num_vertices_left_hemi == num_vertices_right_hemi

wholebrain_coords = np.vstack((left_coords, right_coords))
wholebrain_faces = np.vstack((left_faces, right_faces + num_vertices_left_hemi))