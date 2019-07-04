features_freesurfer = ('freesurfer_thickness',
                       'freesurfer_sulc',
                       'freesurfer_curv',
                       'freesurfer_area',
                       'freesurfer_pial_area',
                       'freesurfer_pial_lgi',
                       'freesurfer_jacobian_white',
                       'freesurfer_volume')
features_fsl = ('gmdensity',)

features_spm_cat_prefixes = { 'spm_cat_gmdensity': 'mwp1',
                              'spm_cat_wmdensity': 'mwp2',
                              }
features_spm_cat = tuple(features_spm_cat_prefixes.keys())

features_cortical = features_freesurfer
features_volumetric = features_fsl + features_spm_cat

base_feature_list = features_cortical + features_volumetric

default_feature_single_edge = ('freesurfer_thickness',)
default_features_multi_edge = ('freesurfer_thickness', 'freesurfer_curv')

default_weight_method = ('manhattan',)

weights_on_original_features = ('diff_medians', 'diff_medians_abs',
                                'diff_means', 'diff_means_abs')

histogram_weights = (
    'chebyshev', 'chebyshev_neg', 'chi_square',
    'correlate', 'correlate_1',
    'cosine', 'cosine_1', 'cosine_2', 'cosine_alt',
    'euclidean', 'fidelity_based',
    'histogram_intersection', 'histogram_intersection_1',
    'jensen_shannon', 'kullback_leibler', 'manhattan', 'minowski',
    'noelle_1', 'noelle_2', 'noelle_3', 'noelle_4', 'noelle_5',
    'relative_bin_deviation', 'relative_deviation')

implemented_weights = histogram_weights + weights_on_original_features

default_minimum_num_bins = 5
default_num_bins = 25
default_trim_percentile = 5

default_atlas = 'fsaverage'  # 'glasser2016'
default_smoothing_param = 10
default_node_size = None

edge_range_predefined = {'freesurfer_thickness': (0.0, 5.0),
                         'freesurfer_curv'     : (-0.3, +0.3),
                         'freesurfer_sulc'     : (-1.5, +1.5),
                         'freesurfer_area'     : (0.0, 1.5)
                         }
default_edge_range = None  # edge_range_predefined[default_feature_single_edge]

default_roi_statistic = 'median'
default_num_procs = 2

# multiedge


# -----------------------------------------------------------------------------------------------
#  atlas and parcellation related
# -----------------------------------------------------------------------------------------------

# keep all the names in lowercase
atlas_list = ['fsaverage', 'glasser2016',
              'yeo2011_fsaverage5', 'yeo2011_fsaverage6', 'yeo2011_fsaverage_highres',
              'cat_aal', 'cat_lpba40', 'cat_ibsr']

default_vbm_atlas = 'cat_aal'

# roi labelled ?? in Glasser parcellation has label 16777215
# fsaverage: label unknown --> 1639705, corpuscallosum --> 3294840
labels_to_ignore_fsaverage_format = [1639705, 3294840]
label_names_to_ignore_fsaverage_format = ['unknown', 'corpuscallosum',
                                          'lh_unknown', 'lh_corpuscallosum',
                                          'rh_unknown', 'rh_corpuscallosum']
ignore_roi_labels = {'glasser2016'              : [16777215, ],
                     'fsaverage'                : labels_to_ignore_fsaverage_format,
                     'yeo2011_fsaverage5'       : labels_to_ignore_fsaverage_format,
                     'yeo2011_fsaverage6'       : labels_to_ignore_fsaverage_format,
                     'yeo2011_fsaverage_highres': labels_to_ignore_fsaverage_format}
ignore_roi_names = {'glasser2016'              : ['??', '???',
                                                  'lh_???', 'rh_???',
                                                  'lh_???', 'rh_???'],
                    'fsaverage'                : label_names_to_ignore_fsaverage_format,
                    'yeo2011_fsaverage5'       : label_names_to_ignore_fsaverage_format,
                    'yeo2011_fsaverage6'       : label_names_to_ignore_fsaverage_format,
                    'yeo2011_fsaverage_highres': label_names_to_ignore_fsaverage_format}

null_roi_index = 0
null_roi_name = 'null_roi_ignore'

# -----------------------------------------------------------------------------------------------

multi_edge_summary_func_default = ('prod', 'median')

if __name__ == '__main__':
    pass
