
__all__ = ['run_workflow', 'extract', 'roiwise_stats_indiv',
           'parcellate', 'freesurfer', 'read_atlas', 'cli_run']

from sys import version_info

if version_info.major==2 and version_info.minor==7:
    import parcellate
    import freesurfer
    from run_workflow import extract, roiwise_stats_indiv, cli_run
    from parcellate import read_atlas
elif version_info.major > 2:
    from graynet import parcellate, freesurfer
    from graynet.run_workflow import extract, roiwise_stats_indiv, cli_run
    from graynet.parcellate import read_atlas
else:
    raise NotImplementedError('hiwenet supports only 2.7.13 or 3+. Upgrade to Python 3+ is recommended.')

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
