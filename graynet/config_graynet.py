
features_freesurfer = ['freesurfer_thickness',
                       'freesurfer_sulc',
                       'freesurfer_curv',
                       'freesurfer_pial_area',
                       'freesurfer_pial_lgi',
                       'freesurfer_jacobian_white',
                       'freesurfer_volume']
features_fsl = ['gmdensity', ]

base_feature_list = features_freesurfer + features_fsl

default_feature_single_edge = ['freesurfer_thickness', ]
default_features_multi_edge = ['freesurfer_thickness', 'freesurfer_curv']

default_weight_method = ('manhattan',)

implemented_weights = [
    'chebyshev', 'chebyshev_neg', 'chi_square',
    'correlate', 'correlate_1',
    'cosine', 'cosine_1', 'cosine_2', 'cosine_alt',
    'euclidean', 'fidelity_based',
    'histogram_intersection', 'histogram_intersection_1',
    'jensen_shannon', 'kullback_leibler', 'manhattan', 'minowski',
    'noelle_1', 'noelle_2', 'noelle_3', 'noelle_4', 'noelle_5',
    'relative_bin_deviation', 'relative_deviation']

default_minimum_num_bins = 5
default_num_bins = 25
default_trim_percentile = 5

default_atlas = 'FSAVERAGE' # 'GLASSER2016'
default_smoothing_param = 10
default_node_size = None

edge_range_predefined = {'freesurfer_thickness': (0, 5),
                         'freesurfer_curv'     : (-0.3, +0.3)}
default_edge_range = None # edge_range_predefined[default_feature_single_edge]

default_roi_statistic = 'median'
default_num_procs = 2

# multiedge

multi_edge_summary_func_default = 'prod' # 'median'

if __name__ == '__main__':
    pass
