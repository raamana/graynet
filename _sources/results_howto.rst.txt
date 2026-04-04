
How to use `graynet` results
----------------------------

graynet 2.0 stores results as run-level Parquet tables plus JSON metadata.
GraphML and CSV are optional export formats, not the canonical on-disk layout.

Loading a run
-------------

.. code-block:: python

    from graynet import load_run

    run = load_run("/path/to/run_dir")
    print(run.metadata.keys())
    print(run.raw_edges)

Filtering raw edges
-------------------

.. code-block:: python

    from graynet import get_edge_values

    edge_rows = get_edge_values(
        run.raw_edges,
        subject_id="sub-001",
        base_feature="freesurfer_thickness",
        weight_method="manhattan",
    )

Converting to an ML matrix
--------------------------

.. code-block:: python

    X = run.raw_edges.to_ndarray(
        ["sub-001", "sub-002", "sub-003"],
        base_feature="freesurfer_thickness",
        weight_method="manhattan",
    )

Reconstructing a NetworkX graph
-------------------------------

.. code-block:: python

    from graynet import export_to_nx

    graph = export_to_nx(edge_rows)

    # centroid coordinates are added to node attributes when available
    print(graph.nodes["lh_bankssts"])

Iterating subject by subject
----------------------------

.. code-block:: python

    for subject_id, subject_edges in run.raw_edges.iter_subjects():
        print(subject_id, subject_edges.to_pandas().shape)

Optional exports
----------------

If you need GraphML or CSV, use the CLI export subcommands documented in
:doc:`usage_cli`.

Canonical tables are a better starting point for:

- machine-learning matrices
- subject-wise filtering
- tabular analysis in Python or R
- reproducible export to GraphML or CSV later
