
Command line interface
-----------------------

The command line interface for graynet (preferred interface, given its target is HPC) is shown below. Check the bottom of this page for examples.

.. argparse::
   :ref: graynet.graynet.__get_parser
   :prog: graynet
   :nodefault:
   :nodefaultconst:

If ``graynet`` is installed without error, you should have ``graynet`` command in your path. Just type ``graynet`` and you should see it display options and usage instructions. For example :

.. code-block:: bash

  $ 10:32:14 SQuark ~ >>  graynet
  usage: graynet [-h] -s SUBJECT_IDS_PATH -i INPUT_DIR [-f FEATURE] [-o OUT_DIR]
                 (-w [WEIGHT_METHOD [WEIGHT_METHOD ...]] | -r [ROI_STATS [ROI_STATS ...]])
                 [-e min max] [-b NUM_BINS] [-a ATLAS] [-n NODE_SIZE]
                 [-p SMOOTHING_PARAM]

  optional arguments:
    -h, --help            show this help message and exit
    -s SUBJECT_IDS_PATH, --subject_ids_path SUBJECT_IDS_PATH
                          Path to file containing list of subject IDs (one per
                          line)
    -i INPUT_DIR, --input_dir INPUT_DIR
                          Path to a folder containing input data. It could ,for
                          example, be a Freesurfer SUBJECTS_DIR, if the chosen
                          feature is from Freesurfer output.
    -f FEATURE, --feature FEATURE
                          Atlas to use to define nodes/ROIs. Default:
                          'freesurfer_thickness'
    -o OUT_DIR, --out_dir OUT_DIR
                          Where to save the extracted features.
    -w [WEIGHT_METHOD [WEIGHT_METHOD ...]], --weight_method [WEIGHT_METHOD [WEIGHT_METHOD ...]]
                          List of methods used to estimate the weight of the
                          edge between the pair of nodes.
    -r [ROI_STATS [ROI_STATS ...]], --roi_stats [ROI_STATS [ROI_STATS ...]]
                          Option to compute summary statistics within each ROI
                          of the chosen parcellation. These statistics (such as
                          the median) can serve as a baseline for network-level
                          values produced by graynet. Options for summary
                          statistics include 'median', 'entropy', 'kurtosis' and
                          any other appropriate summary statistics listed under
                          scipy.stats: https://docs.scipy.org/doc/scipy/referenc
                          e/stats.html#statistical-functions .

  Weight parameters:
    Parameters relevant to histogram edge weight calculations

    -e min max, --edge_range min max
                          The range of edges (two finite values) within which to
                          bin the given values e.g. --edge_range 0.0 5.0 This
                          can be helpful (and important) to ensure
                          correspondence across multiple invocations of graynet
                          (for different subjects), in terms of range across all
                          bins as well as individual bin edges. Default : (0,
                          5), to automatically compute from the given values.
    -b NUM_BINS, --num_bins NUM_BINS
                          Number of bins used to construct the histogram within
                          each ROI or group. Default : 25

  Atlas:
    Parameters describing the atlas, its parcellation and any smoothing of
    features.

    -a ATLAS, --atlas ATLAS
                          Name of the atlas to define parcellation of
                          nodes/ROIs. Default: 'GLASSER2016'
    -n NODE_SIZE, --node_size NODE_SIZE
                          Size of individual node for the atlas parcellation.
                          Default : None
    -p SMOOTHING_PARAM, --smoothing_param SMOOTHING_PARAM
                          Smoothing parameter for feature. Default: FWHM of 10
                          for Freesurfer thickness
  WARNING:root:Too few arguments!


A rough example of usage can be:

.. code-block:: bash

    #!/bin/bash
    #$ -l mf=2G -q queue_name.q -wd /work/project
    cd /work/project
    graynet -s subject_ids.txt -f freesurfer_thickness -i /work/project/freesurfer_reconall -w manhattan eucledian chebyshev -a GLASSER2016 -p 10 -o /work/project/graynet_processing


Note you can specify mulitple weight metrics to save on I/O activity and walltime on HPC.


To automate the the process of generating the above job files and their submission for processing on HPC clusters, I have included a `script <https://github.com/raamana/graynet/blob/master/scripts/generate_hpc_jobs.py>`_ . Please read the instructions inside and run it.

