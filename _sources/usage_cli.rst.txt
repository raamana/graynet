
Command line interface
-----------------------

The command line interface for graynet (preferred interface, given its target is HPC) is shown below. Check the bottom of this page for examples.

.. argparse::
   :ref: graynet.run_workflow.get_parser
   :prog: graynet
   :nodefault:
   :nodefaultconst:


Simple example
---------------

A rough example of usage can be:

.. code-block:: bash

    #!/bin/bash
    #$ -l mf=2G -q queue_name.q -wd /work/project
    cd /work/project
    graynet -s subject_ids.txt -f freesurfer_thickness -i /work/project/freesurfer_reconall \
        -w manhattan eucledian chebyshev -a GLASSER2016 -p 10 -o /work/project/graynet_processing


Note you can specify mulitple weight metrics to save on I/O activity and walltime on HPC.
