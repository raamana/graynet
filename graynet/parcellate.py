"""
Module with handling the parcellation of different cortical atlases.

"""

import os
import os.path
from os.path import isdir
from pathlib import Path

import nibabel as nib
import numpy as np

from graynet import config_graynet as cfg
from graynet.utils import check_atlas, check_atlas_annot_exist, roi_info


def get_atlas_path(atlas_name=None):
    """Validates the atlas name and returns its location"""

    if atlas_name in [None, 'None', '']:
        atlas_name = 'fsaverage'

    atlas_name, _ = check_atlas(atlas_name)

    if atlas_name  in cfg.atlas_list:

        this_dir = Path(__file__).resolve().parent
        if atlas_name in ['glasser2016']:
            atlas_path = this_dir / 'atlases' / 'glasser2016' / \
                         'fsaverage_annot_figshare3498446'
        elif atlas_name in ['fsaverage', 'yeo2011_fsaverage5', 'yeo2011_fsaverage6',
                            'yeo2011_fsaverage_highres']:
            atlas_path = this_dir / 'atlases' /  atlas_name
        elif atlas_name in ['cat_aal', 'cat_lpba40', 'cat_ibsr']:
            # TODO inconsistency: this returns a path to a FILE,
            #   whereas cortical atlases are referred to by their FS folders
            atlas_path = this_dir / 'atlases' /  atlas_name / 'atlas.nii'
        else:
            raise NotImplementedError(
                'Atlas {} is not implemented / unreadable.'.format(atlas_name))

    # cortical atlas in Freesurfer org
    elif os.path.isdir(atlas_name) and check_atlas_annot_exist(atlas_name):
        dir_path_atlas = Path(atlas_name).resolve()
        atlas_path = dir_path_atlas.parent
        atlas_name = dir_path_atlas.name
    else:
        raise NotImplementedError('Invalid choice for atlas!')

    return atlas_path.resolve(), atlas_name


def get_atlas_annot(atlas_name=None):
    """High level wrapper to get all the info just by using a name."""

    atlas_path, atlas_name = get_atlas_path(atlas_name)
    annot = read_atlas_annot(atlas_path)

    return annot, atlas_path


def freesurfer_roi_labels(atlas_name, node_size=None):
    """
    Returns just the vertex-wise indices for grouping the vertices into ROIs.
        Order:  left followed by right.

    """

    annot, _ = get_atlas_annot(atlas_name)
    roi_labels = __combine_annotations(annot, atlas_name)

    if node_size is not None:
        patch_index = load_subdivision_patchwise(atlas_name, node_size).ravel()

        numbered = list()
        for rl, pi in zip(roi_labels, patch_index):
            if rl == cfg.null_roi_name or np.isnan(pi):
                numbered.append(cfg.null_roi_name)
                # as there is no subdivision in a Null ROI
            else:
                numbered.append('{}_p{}'.format(rl, int(pi)))
        # results --> 'lh_precuneus_p9' indicating patch 9 of left precuneus
        #   except when the ROI was labelled cfg.null_roi_name / 'null_roi_ignore'

        roi_labels = np.array(numbered)

    return roi_labels, annot


def __combine_annotations(annot, atlas_name):
    """Combines named labels from two hemispheres, ignoring non-cortex"""

    max_len = 1 + max(max(map(len, annot['lh']['names'] + annot['rh']['names'])),
                      len(cfg.null_roi_name))
    str_dtype = np.dtype('U{}'.format(max_len))

    named_labels = dict()
    for hemi in ['lh', 'rh']:
        named_labels[hemi] = np.empty(annot[hemi]['labels'].shape, str_dtype)
        uniq_labels = np.unique(annot[hemi]['labels'])
        for label in uniq_labels:
            if not ( (label == cfg.null_roi_index) or
                     (label in cfg.ignore_roi_labels[atlas_name])):  # to be ignored
                idx_roi = np.nonzero(annot[hemi]['ctab'][:, 4] == label)[0][0]
                mask_label = annot[hemi]['labels'] == label
                named_labels[hemi][mask_label] = \
                    '{}_{}'.format(hemi, annot[hemi]['names'][idx_roi])

        # setting the non-accessed vertices (part of non-cortex)
        # to a specific label to ignore later
        null_mask = named_labels[hemi] == ''
        named_labels[hemi][null_mask] = cfg.null_roi_name

    # wb stands for whole brain
    wb_named_labels = np.hstack((named_labels['lh'], named_labels['rh']))

    for ignore_label in cfg.ignore_roi_names[atlas_name]:
        wb_named_labels[wb_named_labels == ignore_label] = cfg.null_roi_name

    # # original implementation for glasser2016,
    # #  with labels in different hemi coded differently
    # roi_labels = np.hstack((annot['lh']['labels'], annot['rh']['labels']))

    return wb_named_labels


def read_atlas_annot(atlas_dir, hemi_list=None):
    """ Returns atlas annotations """

    if hemi_list is None:
        hemi_list = ['lh', 'rh']

    annot = dict()
    for hemi in hemi_list:
        annot[hemi] = dict()
        annot_path = atlas_dir / 'label' / '{}.aparc.annot'.format(hemi)
        annot[hemi]['labels'], annot[hemi]['ctab'], \
            annot[hemi]['names'] = nib.freesurfer.io.read_annot(annot_path, orig_ids=True)

        # ensuring names are a plain string
        if isinstance(annot[hemi]['names'][0], np.bytes_):
            annot[hemi]['names'] = \
                [bytestr.decode('UTF-8') for bytestr in annot[hemi]['names']]

    return annot


def read_freesurfer_atlas(atlas_spec, hemi_list=None):
    """Script to read pre-computed parcellations for fsaverage and HCP-MMP-1.0 """

    if hemi_list is None:
        hemi_list = ['lh', 'rh']

    if isdir(atlas_spec):
        atlas_dir = Path(atlas_spec).resolve()
    else:
        atlas_dir, atlas_name = get_atlas_path(atlas_spec)

    coords, faces, info = dict(), dict(), dict()

    for hemi in hemi_list:
        hemi_path = os.path.join(atlas_dir, 'surf', '{}.orig'.format(hemi))
        coords[hemi], faces[hemi], info[hemi] = \
            nib.freesurfer.io.read_geometry(hemi_path, read_metadata=True)

    num_vertices_left_hemi = coords['lh'].shape[0]

    coords['whole'] = np.vstack((coords['lh'], coords['rh']))
    faces['whole'] = np.vstack((faces['lh'], faces['rh'] + num_vertices_left_hemi))

    annot = read_atlas_annot(atlas_dir, hemi_list)

    return coords, faces, annot


def roi_labels_centroids(atlas_name, node_size=None):
    """Returns a list of ROI centroids, for use in visualizations (nodes on a network)"""

    atlas_dir, atlas_name = get_atlas_path(atlas_name)
    coords, faces, annot = read_freesurfer_atlas(atlas_dir)
    vertex_labels, ctx_annot = freesurfer_roi_labels(atlas_name, node_size)
    uniq_rois, roi_size, num_nodes = roi_info(vertex_labels)

    centroids = dict()
    for roi in uniq_rois:
        this_roi_vertices = coords['whole'][vertex_labels == roi, :]
        centroids[roi] = np.median(this_roi_vertices, axis=0)

    return uniq_rois, centroids, vertex_labels


def subdivide_cortex(atlas_dir, hemi_list=None):
    """Subdivides the given cortical parcellation (each label into smaller patches)"""

    raise NotImplementedError('This function has not been implemented yet.')

    # # noinspection PyUnreachableCode
    # if hemi_list is None:
    #     hemi_list = ['lh', 'rh']
    #
    # coords, faces, annot = read_freesurfer_atlas(atlas_dir)
    #
    # labels_to_remove = ['corpuscallosum', 'unknown']
    # null_label = 0
    #
    # def ismember(A, B):
    #     B_unique_sorted, B_idx = np.unique(B, return_index=True)
    #     B_in_A_bool = np.in1d(B_unique_sorted, A, assume_unique=True)
    #
    # cortex_label = dict()
    # for hemi in hemi_list:
    #     cortex_label_path = atlas_dir / 'label' / '{}.cortex.label'.format(hemi)
    #     cortex_label[hemi] = nib.freesurfer.io.read_label(cortex_label_path)
    #
    #     # # cortex_label[hemi] is an index into annot[hemi]['labels']
    #
    #     mask_for_cortex = np.in1d(annot[hemi]['labels'], cortex_label,
    #                               assume_unique=True)


def load_subdivision_patchwise(atlas_name, min_vtx_per_patch=100):
    """
    Loads a precomputed subdivision of the cortex
    Originally generated in Matlab based on k-means clustering (Raamana PhD thesis);
     adaptive subdivision ensuring a specified minimum number of vertices per patch

    more details --> graynet/scripts/export_precomputed_ctx_parc_in_numpy_format.py
    """

    atlas_name = atlas_name.lower()
    if atlas_name not in ('fsaverage',):
        raise ValueError('Only fsaverage is supported at the moment!')

    if min_vtx_per_patch not in cfg.allowed_mvpp:
        raise ValueError('Invalid min_vtx_per_patch. Choose one of {}'
                         ''.format(cfg.allowed_mvpp))

    this_dir = Path(__file__).resolve().parent
    parc_dir = this_dir / 'resources' / 'CorticalSubdivisionIntoPatches' / atlas_name

    file_name = lambda mvp: 'CortexSubdivision_Kmeans_' \
                           'MinVertexCountPerPartition{}.npy' \
                           ''.format(mvp)
    parc = np.load((parc_dir / file_name(min_vtx_per_patch)))

    return parc