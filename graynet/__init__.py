
__all__ = ['graynet', 'extract', 'roiwise_stats_indiv',
           'parcellate', 'freesurfer', 'read_atlas', 'cli_run']

from sys import version_info

if version_info.major==2 and version_info.minor==7:
    import freesurfer
    import parcellate
    from graynet import extract, roiwise_stats_indiv, cli_run
    from parcellate import read_atlas
elif version_info.major > 2:
    from graynet import parcellate, freesurfer
    from graynet.graynet import extract, roiwise_stats_indiv, cli_run
    from graynet.parcellate import read_atlas
else:
    raise NotImplementedError('hiwenet supports only 2.7.13 or 3+. Upgrade to Python 3+ is recommended.')
