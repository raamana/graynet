import networkx as nx
import numpy as np
from graynet.utils import warn_nan
from pathlib import Path

# <<< change the out_dir to the folder contianing the GraphML files here:
out_dir = '/Users/Reddy/dev/graynet/example_data/freesurfer/test_outputs'

out_dir = Path(out_dir).resolve()
for gml_path in out_dir.rglob('*.graphml'):

    graph = nx.read_graphml(gml_path)

    weight_vec = np.array(list(
            nx.get_edge_attributes(graph, 'weight').values()))

    warn_nan(weight_vec)

    # you can control for the format of the file here
    new_path = gml_path.with_suffix('.csv')

    np.savetxt(new_path, weight_vec, fmt='%.5f')