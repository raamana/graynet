
__all__ = ['run_workflow', 'extract', 'roiwise_stats_indiv', 'draw3Dnx',
           'parcellate', 'freesurfer', 'read_freesurfer_atlas', 'cli_run', '__version__']

from ._version import get_versions
__version__ = get_versions()['version']

from sys import version_info
if version_info.major > 2:
    from graynet import utils, parcellate, freesurfer
    from graynet.run_workflow import extract, roiwise_stats_indiv, cli_run
    from graynet.parcellate import read_freesurfer_atlas
    from graynet.vis_network import draw3Dnx
else:
    raise NotImplementedError('graynet supports Python 3 or higher only. Please upgrade.')

del get_versions
del version_info