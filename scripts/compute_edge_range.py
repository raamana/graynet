"""#---------------------------------#---------------------------------#-------------

This script helps estimate the full range of values for a given feature type

"""  #
# ---------------------------------#---------------------------------#-------------

from os.path import join as pjoin

import numpy as np

from graynet.utils import import_features, check_subjects

trim_percentile = 5

in_dir_vbm = '/Users/Reddy/dev/graynet/example_data/volumetric_CAT12'
sub_list, _, _, _ = check_subjects(pjoin(in_dir_vbm, 'sub_id_list.txt'))

base_feature = 'spm_cat_gmdensity'


features = import_features(in_dir_vbm, sub_list, base_feature)
print('\nThe min max per subject (after trimming {}% outliers):'.format(
        trim_percentile))

min_all =  np.Inf
max_all = -np.Inf
for sid, data in features.items():
    vec = np.array(data).flatten()
    min_sub = np.percentile(vec, trim_percentile)
    max_sub = np.percentile(vec, 100 - trim_percentile)

    print('{} {} {}'.format(sid, min_sub, max_sub))
    if min_sub < min_all:
        min_all = min_sub

    if max_sub > max_all:
        max_all = max_sub

print('\nAfter computing subject-wise min and max after trimming {}% outliers,\n'
      '\tthe dataset-wide min and max are: {} and {}'
      ''.format(trim_percentile,min_all, max_all))
