
Command line interface
-----------------------

graynet 2.0 exposes subcommands for the main run modes and for optional exports.
Each successful command writes one run directory containing canonical outputs.

.. argparse::
   :ref: graynet.cli.build_parser
   :prog: graynet
   :nodefault:
   :nodefaultconst:


Single-feature edge extraction
------------------------------

.. code-block:: bash

    graynet edges \
      -i /work/project/freesurfer_reconall \
      -s subject_ids.txt \
      -f freesurfer_thickness \
      -w manhattan euclidean \
      -a fsaverage \
      -p 10 \
      -o /work/project/graynet_runs \
      -c 2

Multi-edge extraction
---------------------

.. code-block:: bash

    graynet multiedge \
      -i /work/project/freesurfer_reconall \
      -s subject_ids.txt \
      -f freesurfer_thickness freesurfer_curv \
      -w manhattan \
      -t median prod \
      -a fsaverage \
      -p 10 \
      -o /work/project/graynet_runs

ROI-wise statistics
-------------------

.. code-block:: bash

    graynet roi-stats \
      -i /work/project/freesurfer_reconall \
      -s subject_ids.txt \
      -f freesurfer_thickness \
      -r median mean \
      -a fsaverage \
      -p 10 \
      -o /work/project/graynet_runs

Exporting GraphML or CSV
------------------------

Canonical outputs are Parquet plus JSON metadata. Export GraphML or CSV only
when you need those downstream formats.

.. code-block:: bash

    graynet export graphml --run-dir /work/project/graynet_runs/edges__atlas-fsaverage__features-freesurfer_thickness__weights-manhattan__smth-10__node-none

.. code-block:: bash

    graynet export csv --run-dir /work/project/graynet_runs/edges__atlas-fsaverage__features-freesurfer_thickness__weights-manhattan__smth-10__node-none
