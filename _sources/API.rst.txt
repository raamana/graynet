--------------
API
--------------

A tutorial-like presentation is available at :doc:`usage_api`.

The core public API is:

- ``graynet.extract`` for single-feature edge runs
- ``graynet.extract_multiedge`` for multi-feature runs with summary edges
- ``graynet.roiwise_stats_indiv`` for ROI-wise summary statistics
- ``graynet.load_run`` and related helpers for reading canonical run outputs

graynet 2.0 writes canonical run-level outputs in Parquet plus JSON metadata.
GraphML and CSV are now optional exports derived from those canonical tables.

API Reference
--------------

.. automodule:: graynet.api
   :members: extract, roiwise_stats_indiv, extract_multiedge
   :undoc-members:
   :show-inheritance:

.. automodule:: graynet.results
   :members: RunData, EdgeData, RoiStatsData, load_run, get_edge_values, export_to_nx
   :undoc-members:
   :show-inheritance:
