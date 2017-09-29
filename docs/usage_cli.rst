
Command line interface
-----------------------

The command line interface for graynet (preferred interface, given its target is HPC) is shown below. Check the bottom of this page for examples.

.. argparse::
   :ref: graynet.run_workflow.__get_parser
   :prog: graynet
   :nodefault:
   :nodefaultconst:


A rough example of usage can be:

.. code-block:: bash

    #!/bin/bash
    #$ -l mf=2G -q queue_name.q -wd /work/project
    cd /work/project
    graynet -s subject_ids.txt -f freesurfer_thickness -i /work/project/freesurfer_reconall -w manhattan eucledian chebyshev -a GLASSER2016 -p 10 -o /work/project/graynet_processing


Note you can specify mulitple weight metrics to save on I/O activity and walltime on HPC.


To automate the the process of generating the above job files and their submission for processing on HPC clusters, I have included a `script <https://github.com/raamana/graynet/blob/master/scripts/generate_hpc_jobs.py>`_ . Please read the instructions inside and run it.

