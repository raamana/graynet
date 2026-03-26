
Command line interface
-----------------------

The command line interface for graynet (preferred interface, given its target is HPC) is shown below. Check the bottom of this page for examples.

.. argparse::
   :ref: graynet.cli.build_parser
   :prog: graynet
   :nodefault:
   :nodefaultconst:


A rough example of usage can be:

.. code-block:: bash

    #!/bin/bash
    #$ -l mf=2G -q queue_name.q -wd /work/project
    cd /work/project
    graynet edges -s subject_ids.txt -f freesurfer_thickness -i /work/project/freesurfer_reconall -w manhattan euclidean chebyshev -a glasser2016 -p 10 -o /work/project/graynet_processing


Note you can specify mulitple weight metrics to save on I/O activity and walltime on HPC.
