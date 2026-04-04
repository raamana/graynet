------------------
Examples using API
------------------

graynet 2.0 keeps the classic entrypoints for computing outputs, but the result
layout is now run-level and canonical:

- ``run_metadata.json``
- ``edges_raw.parquet``
- ``edges_summary.parquet`` when summary statistics are requested
- ``roi_stats.parquet`` when ROI statistics are requested

Single-feature edge extraction
------------------------------

.. code-block:: python

    from pathlib import Path
    from graynet import extract

    input_dir = Path("/project/freesurfer_reconall")
    subjects = ["sub-001", "sub-002"]

    run_dir = extract(
        subjects,
        input_dir,
        base_feature="freesurfer_thickness",
        weight_method_list=["manhattan", "cosine"],
        atlas="fsaverage",
        smoothing_param=10,
        out_dir=Path("/project/graynet_runs"),
        return_results=False,
        num_procs=2,
    )

    print(run_dir)

Reading run outputs
-------------------

.. code-block:: python

    from graynet import load_run, get_edge_values, export_to_nx

    run = load_run(run_dir)
    print(run.metadata["subject_ids"])
    print(run.metadata["node_labels"][:5])

    edge_rows = get_edge_values(
        run.raw_edges,
        subject_id="sub-001",
        base_feature="freesurfer_thickness",
        weight_method="manhattan",
    )

    X = run.raw_edges.to_ndarray(
        ["sub-001", "sub-002"],
        base_feature="freesurfer_thickness",
        weight_method="manhattan",
    )

    graph = export_to_nx(edge_rows)

Multi-edge extraction
---------------------

.. code-block:: python

    from graynet import extract_multiedge

    run_dir = extract_multiedge(
        subjects,
        input_dir,
        base_feature_list=["freesurfer_thickness", "freesurfer_curv"],
        weight_method_list=["manhattan"],
        summary_stats=["median", "prod"],
        atlas="fsaverage",
        smoothing_param=10,
        out_dir=Path("/project/graynet_runs"),
        return_results=False,
    )

ROI-wise statistics
-------------------

.. code-block:: python

    import numpy as np
    from graynet import roiwise_stats_indiv

    run_dir = roiwise_stats_indiv(
        subjects,
        input_dir,
        base_feature="freesurfer_thickness",
        chosen_roi_stats=["median", np.nanmean],
        atlas="fsaverage",
        smoothing_param=10,
        out_dir=Path("/project/graynet_runs"),
        return_results=False,
    )




