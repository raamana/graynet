import os
import shlex
import sys
from os.path import abspath, dirname, exists as pexists, join as pjoin, realpath
from sys import version_info

sys.dont_write_bytecode = True

if __name__ == '__main__' and __package__ is None:
    parent_dir = dirname(dirname(abspath(__file__)))
    pkg_dir = dirname(parent_dir)
    sys.path.append(parent_dir)
    sys.path.append(pkg_dir)

if version_info.major > 2:
    from graynet.run_workflow import cli_run as CLI
else:
    raise NotImplementedError('graynet requires Python 3+.')

test_dir = dirname(os.path.realpath(__file__))
base_dir = realpath(pjoin(test_dir, '..', '..', 'example_data'))

subject_id_list = ['subject12345', ]

in_dir_vbm = pjoin(base_dir, 'volumetric_CAT12')
base_feature = 'spm_cat_gmdensity'
# atlas = 'CAT_LPBA40' # 'CAT_IBSR' # 'CAT_AAL'  #

cur_dir = os.path.dirname(abspath(__file__))
atlas_dir = abspath(pjoin(cur_dir, '..', 'atlases'))
atlas = pjoin(atlas_dir, 'cat_aal', 'aal.nii')

sub_list = pjoin(in_dir_vbm, 'sub_id_list_test.txt')

out_dir = pjoin(in_dir_vbm, 'vol_CAT12_graynet')
if not pexists(out_dir):
    os.mkdir(out_dir)


def test_CLI_weight():
    """ensures the CLI works. """

    sys.argv = shlex.split('graynet -i {} -s {} -f {} -w manhattan -o {} -a {} -c 1'
                           ''.format(in_dir_vbm, sub_list, base_feature, out_dir,
                                     atlas))

    CLI()


def test_CLI_diff_atlases_by_path():
    """testing for different atlases"""

    for atlas_name in ('aal', 'lpba40', 'ibsr'):
        atlas_path = pjoin(atlas_dir, 'cat_{}'.format(atlas_name),
                           '{}.nii'.format(atlas_name))
        sys.argv = shlex.split('graynet -i {} -s {} -f {} -w manhattan '
                               '-o {} -a {} -c 1'
                               ''.format(in_dir_vbm, sub_list, base_feature,
                                         out_dir, atlas_path))
        CLI()
