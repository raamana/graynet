
from graynet import parcellate
import networkx as nx
import plotly
from os.path import join as pjoin, isdir, realpath, dirname, basename
# import plotly.plotly as py
# plotly.offline.init_notebook_mode(connected=False)

from graynet.vis_network import draw3Dnx

coords, faces, annot = parcellate.read_freesurfer_atlas('fsaverage')

this_dir = dirname(realpath(__file__))
test_out_dir = realpath(pjoin('..', 'example_data', 'freesurfer', 'test_outputs'))
path_graphml = pjoin('subject12345',
                     'freesurfer_thickness_fsaverage_smoothing10_sizeNone_'
                     'edgeweight_diff_medians_abs_graynet.graphml')
brain = nx.read_graphml(path_graphml)

# random=nx.random_geometric_graph(nx.nodes(brain), 0.05)
# # pos=nx.get_node_attributes(random,'pos')
#
# node_attrs = ['x', 'y', 'z']
# for attr in node_attrs:
#     nx.set_node_attributes(random, nx.get_node_attributes(brain, attr), attr)
#
# edge_attrs = ['weight', ]
# for attr in edge_attrs:
#     nx.set_edge_attributes(random, nx.get_edge_attributes(brain, attr), attr)


fig = draw3Dnx(brain, perc_threshold=95)
plotly.offline.plot(fig, filename='networkx.html')