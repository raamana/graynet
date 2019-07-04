"""

Module for all things voxelwise and volumetric

"""

import traceback
import warnings
import hiwenet
import numpy as np
import nibabel
import networkx as nx
import sys

from graynet.utils import import_features, warn_nan, stamp_expt_weight, \
    save_per_subject_graph, save, roi_info, mask_background_roi
from graynet import config_graynet as cfg
from graynet.parcellate import get_atlas_path


def extract_per_subject_volumetric(input_dir, base_feature, roi_labels,
                                   centroids, weight_method_list, atlas,
                                   smoothing_param, node_size, num_bins,
                                   edge_range, out_dir, return_results,
                                   pretty_print_options, subject=None):
    # purposefully leaving subject parameter last to enable partial function creation
    """
    Extracts give set of weights for one subject.

    Parameters
    ----------
    subject
    input_dir
    base_feature
    roi_labels
    weight_method_list
    atlas
    smoothing_param
    node_size
    num_bins
    edge_range
    out_dir
    return_results
    pretty_print_options

    Returns
    -------

    """

    if subject is None:
        return

    try:
        features = import_features(input_dir,
                                   [subject, ],
                                   base_feature,
                                   fwhm=smoothing_param,
                                   atlas=atlas)
    except:
        traceback.print_exc()
        warnings.warn('Unable to read {} features for {}\n Skipping it.'.format(
                base_feature, subject), UserWarning)
        return

    data, rois = mask_background_roi(features[subject], roi_labels,
                                     cfg.null_roi_index) #note:its not null_roi_name!

    max_id_width, nd_id, num_weights, max_wtname_width, nd_wm = pretty_print_options

    if return_results:
        edge_weights_all = dict()
    else:
        edge_weights_all = None

    for ww, weight_method in enumerate(weight_method_list):
        # unique stamp for each subject and weight
        expt_id = stamp_expt_weight(base_feature, atlas, smoothing_param, node_size,
                                    weight_method)
        sys.stdout.write(
            '\nProcessing id {:{id_width}} -- weight {:{wtname_width}} '
            '({:{nd_wm}}/{:{nd_wm}})'
            ' :'.format(subject, weight_method, ww + 1, num_weights,
                        nd_id=nd_id, nd_wm=nd_wm,
                        id_width=max_id_width, wtname_width=max_wtname_width))

        # actual computation of pair-wise features
        try:
            graph = hiwenet.extract(data,
                                    rois,
                                    weight_method=weight_method,
                                    num_bins=num_bins,
                                    edge_range=edge_range,
                                    return_networkx_graph=True)

            # retrieving edge weights
            weight_vec = np.array(list(nx.get_edge_attributes(graph, 'weight').values()))
            warn_nan(weight_vec)
            # weight_vec = get_triu_handle_inf_nan(edge_weights)

            # adding position info to nodes (for visualization later)
            for roi in centroids:
                graph.node[roi]['x'] = float(centroids[roi][0])
                graph.node[roi]['y'] = float(centroids[roi][1])
                graph.node[roi]['z'] = float(centroids[roi][2])

            if return_results:
                edge_weights_all[(weight_method, subject)] = weight_vec

            # saving to disk
            try:
                save(weight_vec, out_dir, subject, expt_id)
                save_per_subject_graph(graph, out_dir, subject, expt_id)
            except:
                raise IOError('Unable to save the network or vectorized weights '
                               'to:\n{}'.format(out_dir))

        except (RuntimeError, RuntimeWarning) as runexc:
            print(runexc)
        except KeyboardInterrupt:
            print('Exiting on keyborad interrupt! \n'
                  'Abandoning the remaining processing for {} weights:\n'
                  '{}.'.format(num_weights - ww, weight_method_list[ww:]))
            sys.exit(1)
        except:
            print('Unable to extract {} features for {}'.format(weight_method, subject))
            traceback.print_exc()

        sys.stdout.write('Done.')

    return edge_weights_all


def volumetric_roi_info(atlas_name):
    """Returns a list of unique ROIs, their labels and centroids"""

    atlas_path, atlas_name = get_atlas_path(atlas_name)
    atlas_labels = nibabel.load(atlas_path).get_data()

    # TODO names for ROIs are not read and used!

    uniq_rois, roi_size, num_nodes = roi_info(atlas_labels, freesurfer_annot=False)

    centroids = dict()
    for roi in uniq_rois:
        centroids[roi] = np.median(np.nonzero(atlas_labels==roi), axis=1)

    return uniq_rois, centroids, atlas_labels
