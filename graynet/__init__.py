__all__ = ['extract', 'roiwise_stats_indiv', 'extract_multiedge', 'load_run',
           'get_edge_values', 'export_to_nx', 'parcellate',
           'freesurfer', 'read_freesurfer_atlas', '__version__']

from importlib.metadata import PackageNotFoundError, version

from graynet import parcellate, freesurfer
from graynet.api import extract, roiwise_stats_indiv, extract_multiedge
from graynet.parcellate import read_freesurfer_atlas
from graynet.results import load_run, get_edge_values, export_to_nx

try:
    __version__ = version('graynet')
except PackageNotFoundError:
    __version__ = '2.0.0'
