import nilearn
import os
import sys
import stat
from os.path import join as pjoin, exists as pexists, realpath
import numpy as np
from matplotlib import pyplot as plt
import pickle
import traceback

from graynet import parcellate
from graynet.run_workflow import roi_info

from nilearn.plotting import plot_connectome

#---------------------------------
# CHANGE THESE FOR YOUR PROCESSING
#---------------------------------

base_dir = '/u1/work/hpc3194'
dataset_name = 'PPMI' # '4RTNI' # 'PPMI' #

# list_of_datasets = [ '4RTNI', 'PPMI', 'ADNI' ]
# list_of_subject_lists = ['graynet.compute.list']*3

proc_dir = pjoin(base_dir, dataset_name, 'processed')
freesurfer_dir = pjoin(proc_dir, 'freesurfer')
target_list_dir = pjoin(proc_dir, 'target_lists')

subject_id_list = pjoin(target_list_dir, 'graynet.compute.list')
id_list = np.atleast_1d(np.genfromtxt(subject_id_list, dtype=str).astype(str))

id_list = id_list[:49]
num_subjects = len(id_list)

base_feature =  'freesurfer_curv' # 'freesurfer_thickness' # 'freesurfer_thickness'
atlas = 'GLASSER2016' # 'FSAVERAGE' # 'GLASSER2016' #
fwhm = 10
this_dir = os.path.dirname(os.path.realpath(__file__))
# atlas_path = os.path.realpath(pjoin(this_dir, '..', 'graynet',
#                                     'atlases', 'glasser2016', 'fsaverage_annot_figshare3498446'))

print('using fsaverage coords for visualization, as Glasser2016 is resampled to it anyways')
atlas_path = os.path.realpath(pjoin(this_dir, '..', 'graynet', 'atlases', 'fsaverage'))
coords, faces, annot = parcellate.read_atlas(atlas_path)

roi_labels, annot = parcellate.freesurfer_roi_labels(atlas)
uniq_rois, roi_size, num_nodes = roi_info(roi_labels)

mean_coords = np.zeros([len(uniq_rois), 3], np.float16)
for rr, roi in enumerate(uniq_rois):
    index_roi = roi_labels == roi
    mean_coords[rr,:] = np.mean(coords['whole'][index_roi], axis=0)

out_dir = pjoin(proc_dir, 'graynet', '{}_{}_fwhm{}'.format(base_feature, atlas, fwhm))

# exclude kullback_leibler as its got too many problems
histogram_dist = np.array([
    'chebyshev', 'chebyshev_neg', 'chi_square',
    'correlate', 'correlate_1',
    'cosine', 'cosine_1', 'cosine_2', 'cosine_alt',
    'euclidean', 'fidelity_based',
    'histogram_intersection', 'histogram_intersection_1',
    'jensen_shannon', 'kullback_leibler', 'manhattan', 'minowski',
    'noelle_1', 'noelle_2', 'noelle_3', 'noelle_4', 'noelle_5',
    'relative_bin_deviation', 'relative_deviation'])
weight = 'chi_square'
num_weights = len(histogram_dist)

def get_adjacency_matrix(out_dir, sid, expt_id):
    "Returns the adjacency matrix"

    vec_path = pjoin(out_dir, sid, '{}_graynet.csv'.format(expt_id))
    edge_vec = np.genfromtxt(vec_path)

    matrix_size = np.int64( (1.0 + np.sqrt(1.0+8.0*len(edge_vec)))/2.0 )
    edge_mat = np.zeros([matrix_size, matrix_size])

    # making this symmetric as required by nilearn's plot_connectome (stupid)
    # upper tri; diag +1; # lower tri; diag -1
    upper_tri = np.triu_indices_from(edge_mat, +1)
    lower_tri = np.tril_indices_from(edge_mat, -1)
    edge_mat[upper_tri] = edge_vec
    edge_mat[lower_tri] = edge_mat.T[lower_tri]

    return edge_mat

def stamp_experiment(base_feature, atlas, smoothing_param, node_size, weight_method):
    "Constructs a string to uniquely identify a given feature extraction method."

    # expt_id = 'feature_{}_atlas_{}_smoothing_{}_size_{}'.format(base_feature, atlas, smoothing_param, node_size)
    expt_id = '{}_{}_smoothing{}_size{}_edgeweight_{}'.format(base_feature, atlas, smoothing_param, node_size, weight_method)

    return expt_id


node_subset = np.hstack( [ np.array(range(20)), np.array(range(359,240,-1)) ] )

saved_path = pjoin(out_dir,'adj_matrices_all_subjects_weights.pkl')
try:
    if pexists(saved_path):
        print('reloading data')
        with open(saved_path,'rb') as sf:
            data_list = pickle.load(sf)
            adj_mat_all = data_list[0]
    else:
        raise IOError('saved data doesnt exist')
except:
    print('re-reading data')
    adj_mat_all = np.zeros([num_subjects, num_weights, num_nodes, num_nodes])
    adj_mat_avg = np.zeros([num_nodes, num_nodes])
    for ss, sid in enumerate(id_list):
        for ww, weight in enumerate(histogram_dist):
            expt_id = stamp_experiment(base_feature, atlas, fwhm, None, weight)
            try:
                adj_mat = get_adjacency_matrix(out_dir, sid, expt_id)
                # replacing nan's with zeros for now
                adj_mat[np.isnan(adj_mat)] = 0
                adj_mat_all[ss, ww, :, :] = adj_mat
            except:
                print('--data for pair-wise dist {} for subject {} not available'.format(weight, sid))

    with open(saved_path,'wb') as sf:
        pickle.dump([adj_mat_all, ], sf)

vis_out_dir = pjoin(out_dir, 'vis')
if not pexists(vis_out_dir):
    os.mkdir(vis_out_dir)

num_rows = 6
num_cols = 4
figsize = [ 14, 18]

def save_fig_get_new(cur_fig, weight, fig_count, img):

    out_fig_path = pjoin(vis_out_dir, '{}_fig{}_{}subjects.pdf'.format(weight, fig_count, num_rows * num_cols))

    cax = cur_fig.add_axes([0.93, 0.2, 0.02, 0.6])
    cur_fig.colorbar(img, cax=cax)
    plt.suptitle(weight, fontsize=20)
    plt.savefig(out_fig_path, dpi=200)
    plt.close()

    new_fig, new_axes = plt.subplots(num_rows, num_cols, figsize=figsize)

    fig_count = fig_count + 1

    return new_fig, new_axes, fig_count


for ww, weight in enumerate(histogram_dist):

    fig_count = 1
    fig, ax = plt.subplots(num_rows, num_cols, figsize=figsize)
    print('visualizing {} pair-wise dist'.format(weight))

    all_subjects_data = np.squeeze(adj_mat_all[:, ww, :, :]).flatten()
    clim_max = np.percentile(all_subjects_data, 99.5)
    clim_min = np.percentile(all_subjects_data, 0.1)

    for ss in range(num_subjects):
        sid = id_list[ss]
        rem = np.mod(ss, num_rows*num_cols)

        try:
            adj_mat = np.squeeze(adj_mat_all[ss, ww, :, :])
            if ss > 0 and rem == 0:
                fig, ax, fig_count = save_fig_get_new(fig, weight, fig_count, adj_mat)

            ax = plt.subplot(num_rows, num_cols, rem + 1)

            img = ax.imshow(adj_mat, vmin=clim_min, vmax=clim_max)

            plt.title('{}'.format(sid))
            ax.xaxis.set_ticklabels([])

            # plot_connectome(adj_mat[:50,:50], mean_coords[:50,:])
            # edge_thr = np.percentile(adj_mat.flatten(), 99)
            # plot_connectome(adj_mat, mean_coords,edge_threshold='99.99%', node_size=5)
        except:
            print('unable to visualize data for pair-wise dist {} for subject {}'.format(weight, sid))
            traceback.print_exc()

    # plt.show(block=False)
    fig, ax, fig_count = save_fig_get_new(fig, weight, fig_count, img)

print('')

