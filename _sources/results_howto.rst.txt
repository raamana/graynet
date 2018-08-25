
How to use `graynet` results
----------------------------


For each combination of the parameters chosen, such as edge metric, atlas etc, ``graynet`` produces one network for each subject. The output format is in ``GraphML`` format, which can be easily read with ``networkx``:

.. code-block:: python

    import networkx as nx
    graph = nx.read_graphml(path_to_graphml_file)


More info can be found `here <https://networkx.github.io/documentation/stable/reference/readwrite/generated/networkx.readwrite.graphml.read_graphml.html>`_.


The graph inside the ``graphML`` file is essentially a pair-wise distance matrix (measured by the metric chosen). There are many ways you can use it - easiest among them is to extract the upper triangular part of the connectivity matrix (as it is symmetric) and treat it a single-subject feature vector for that subject.

These subject-wise feature vectors can be used in many applications, including in the study of brain-behaviour relationships and as a biomarker candidate (e.g. see `this study on ADNI and ABIDE <https://www.biorxiv.org/content/early/2017/07/31/170381>`_). If you are interested in evaluating their predictive utility (out of sample prediction via cross-validation), it's quite simple via `neuropredict <https://github.com/raamana/neuropredict>`_.

We plan to include additional scripts and convenience methods into ``graynet`` to gather the results into readily usable data structures such as `pyradigm <https://github.com/raamana/pyradigm>`_ (or CSV files) for further analysis. Stay tuned!
