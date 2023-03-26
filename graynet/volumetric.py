"""

Module for all things voxelwise and volumetric

"""

import sys
import traceback
import warnings

import hiwenet
import networkx as nx
import nibabel
import numpy as np

from graynet import config_graynet as cfg
from graynet.parcellate import get_atlas_path
from graynet.utils import (import_features, is_image, is_image_3D,
                           mask_background_roi, roi_info, save,
                           save_per_subject_graph, stamp_expt_weight,
                           warn_nan)


def extract_per_subject_volumetric(input_dir, base_feature, roi_labels, centroids,
                                   weight_method_list, atlas_spec, atlas_name,
                                   smoothing_param, node_size, num_bins, edge_range,
                                   out_dir, return_results, pretty_print_options,
                                   subject=None):
    # purposefully leaving subject parameter last to enable partial function creation
    """Extracts a given set of weights for one subject."""

    # identical to extract_per_subject_cortical, except for the use of null roi index
    #   volumetric uses cfg.null_roi_index, different from cfg.null_roi_name
    if subject is None:
        return

    try:
        features = import_features(input_dir,
                                   [subject, ],
                                   base_feature,
                                   fwhm=smoothing_param,
                                   atlas=atlas_spec)
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
        # TODO need a way to identify the atlas by a shot string
        #   when it is originally supplied as an image or a path to an image
        expt_id = stamp_expt_weight(base_feature, atlas_name, smoothing_param,
                                    node_size, weight_method)
        sys.stdout.write(
            '\nProcessing {sid:{id_width}} -- weight {wm:{wtname_width}} '
            '({wc:{nd_wm}}/{nw:{nd_wm}}) :\n'
            ''.format(sid=subject, wm=weight_method, wc=ww + 1, nw=num_weights,
                      nd_id=nd_id, nd_wm=nd_wm, id_width=max_id_width,
                      wtname_width=max_wtname_width))

        # actual computation of pair-wise features
        try:
            graph = hiwenet.extract(data,
                                    rois,
                                    weight_method=weight_method,
                                    num_bins=num_bins,
                                    edge_range=edge_range,
                                    return_networkx_graph=True)

            # retrieving edge weights
            weight_vec = np.array(list(
                    nx.get_edge_attributes(graph, 'weight').values()))
            warn_nan(weight_vec)

            # adding position info to nodes (for visualization later)
            for roi in centroids:
                graph.nodes[roi]['x'] = float(centroids[roi][0])
                graph.nodes[roi]['y'] = float(centroids[roi][1])
                graph.nodes[roi]['z'] = float(centroids[roi][2])

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
            print('Exiting on keyboard interrupt! \n'
                  'Abandoning the remaining processing for {} weights:\n'
                  '{}.'.format(num_weights - ww, weight_method_list[ww:]))
            sys.exit(1)
        except:
            print('Unable to extract {} features for {}'
                  ''.format(weight_method, subject))
            traceback.print_exc()

        sys.stdout.write('Done.')

    return edge_weights_all


def volumetric_roi_info(atlas_spec):
    """Returns a list of unique ROIs, their labels and centroids"""

    if is_image(atlas_spec) and is_image_3D(atlas_spec):
        if atlas_spec.__class__ in nibabel.all_image_classes:
            atlas_labels = atlas_spec.get_fdata()
        else:
            atlas_labels = np.array(atlas_spec)
    elif isinstance(atlas_spec, str):
        atlas_path, atlas_name = get_atlas_path(atlas_spec)
        atlas_labels = nibabel.load(atlas_path).get_fdata()
    else:
        raise ValueError('Unrecognized atlas specification!'
                         'Must be a predefined name, or'
                         'a preloaded image!')

    # TODO names for ROIs are not read and used!

    uniq_rois, roi_size, num_nodes = roi_info(atlas_labels, freesurfer_annot=False)

    centroids = dict()
    for roi in uniq_rois:
        centroids[roi] = np.median(np.nonzero(atlas_labels==roi), axis=1)

    return uniq_rois, centroids, atlas_labels
