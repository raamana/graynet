
from sys import version_info

if version_info.major==2 and version_info.minor==7 and version_info.micro==13:
    import freesurfer
    import parcellate
    from graynet import extract
    from parcellate import read_atlas
elif version_info.major > 2:
    from graynet import parcellate, freesurfer
    from graynet.graynet import extract
    from graynet.parcellate import read_atlas
else:
    raise NotImplementedError('hiwenet supports only 2.7.13 or 3+. Upgrade to Python 3+ is recommended.')

__all__ = ['graynet', 'extract', 'parcellate', 'freesurfer', 'read_atlas']