graynet
===================================

.. image:: http://joss.theoj.org/papers/10.21105/joss.00924/status.svg
   :target: https://doi.org/10.21105/joss.00924

``graynet`` computes single-subject morphometric networks and ROI statistics
from cortical Freesurfer outputs and supported volumetric feature maps.

graynet 2.0 writes canonical run-level outputs as Parquet tables plus JSON
metadata:

- ``run_metadata.json``
- ``edges_raw.parquet``
- ``edges_summary.parquet`` when summary edges are requested
- ``roi_stats.parquet`` when ROI statistics are requested

GraphML and CSV remain available as optional exports.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage_cli
   usage_api
   results_howto
   cortical
   volumetric
   API
   citation


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
