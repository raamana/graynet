Volumetric graynet
------------------

For the volumetric stream of processing, ``graynet`` relies on outputs from `CAT12 toolbox <http://www.neuro.uni-jena.de/cat/>`_ within the SPM ecosystem. *We plan to support various other VBM toolboxes soon* - make your suggestions by `opening an issue <https://github.com/raamana/graynet/issues/new>`_.

The following steps should help you get started and going quickly:

    - Install `CAT12 toolbox <http://www.neuro.uni-jena.de/cat/>`_ and read the `manual <http://www.neuro.uni-jena.de/cat12/CAT12-Manual.pdf>`_
    - Ensure you produce the outputs in *normalized* atlas space, by selecting the output option to be either ``Modulated normalized`` or ``normalized`` depending on your analyses plan. Correct processing would produce images of size 121 x 145 x 121, which must be the same as the atlases included in the CAT12 toolbox
    - Run the CAT12 processing and ensure it is complete
        - Perform the QA as suggested in their `manual <http://www.neuro.uni-jena.de/cat12/CAT12-Manual.pdf>`_ e.g. "SECOND MODULE: DISPLAY ONE SLICE FOR ALL IMAGES".
        - For a more comprehensive QC of segmentation accuracy, you can employ `visualqc <https://raamana.github.io/visualqc/readme.html>`_
    - If everything went well, you would have an ``mri`` folder created in the output folder you chose. This folder would have several files, for each subject you processed e.g. ``mwp1SUBJECTID.nii`` contains the voxel-wise gray matter densities for SUBJECTID in the modulated (prefix ``m`` ) normalized (``w`` ) atlas space. Similarly ``mwp2*.nii`` files would contain voxel-wise white matter densities.


If graynet is installed without error, you should have ``graynet`` command in your path. Just type ``graynet`` and you should see it display options and usage instructions, as presented in :doc:`usage_cli` page.

Suppose say

  - you already have CAT12 processing done at ``/work/project/cat12_vbm``
  - for 5 subjects with identifiers are stored in subject_ids.txt (one per line)
  - you want to analyze network features based on GM density values (identified in ``graynet`` by ``spm_cat_gmdensity`` ), which have been resampled to the space of AAL atlas (identified in ``graynet`` by ``cat_aal`` ) provided in CAT12 (dimensions: 121 x 145 x 121, with voxel sizes: 1.5000, 1.5000, 1.5000)
  - use the ``manhattan`` method to compute histogram distance between two cortical ROIs (assuming all subjects resampled to same atlas)
  - and store the networks extracted in ``/work/project/cat12_vbm/graynet_volumetric``

You could achieve it all with a single command:

.. code-block:: bash

    cd /work/project
    graynet -s subject_ids.txt -f spm_cat_gmdensity -i /work/project/cat12_vbm -w manhattan -a cat_aal -o /work/project/cat12_vbm/graynet_volumetric

That's it! By the time you can get your coffee or stretch your legs, you should have graynet processing done.

Suppose, you prefer to analyze ROIs as defined by other atlases, graynet also supports the IBSR and LPBA40 as supplied by the CAT12 toolbox. If you are also curious to try multiple weight metrics -  ``chebyshev`` and ``cosine``, you can achieve it via the following command:

.. code-block:: bash

    cd /work/project
    graynet -s subject_ids.txt -f spm_cat_gmdensity -i /work/project/cat12_vbm -w chebyshev cosine manhattan -a cat_LPBA40 -o /work/project/cat12_vbm/graynet_volumetric


You could also study other features by simply adding more features to the ``-f``, such as ``spm_cat_wmdensity`` (white matter density).

However, if you have 100s of subjects or your atlas has a large number of ROIs (like the AAL atlas which has 120), this computation can take a long time. It is advisable to move the processing to a larger server or computer cluster, by delegating few subjects to a single job. You could leverage this `script <https://github.com/raamana/graynet/blob/master/scripts/generate_hpc_jobs.py>`_ to process any large dataset at once (produce jobs for processing on a cluster).

    - If you have access to a computer cluster, you can submit the generated jobs to the cluster.
    - If you don't have access to a cluster, you could still use the generated job/script files to run them locally - few at a time.

After processing is done (on the cluster or locally), and assuming no other errors, you should see a CSV file for each subject, with ``n*(n-1)/2`` values corresponding to the upper-triangular part of the connectivity matrix, where n is the number of ROIs in a atlas you chose (68 for fsaverage or 122 for AAL). In addition to easy to use CSV files, graynet also produces rich graph structure (in the ``graphml`` format easily usable via the ``networkx`` package) that identifies different nodes in these subject-wise network and the edges for each node, along with their edge weights.

The output folder will be graynet within ``proc_dir`` you choose in the script - determined by this statement in the above script:

.. code-block:: python

    out_dir = pjoin(proc_dir, 'graynet', '{}_{}_range{}_{}_nbins{}'.format(base_feature, atlas, edge_range[0], edge_range[1], num_bins))


which after expansion looks something like based on your choices: ``/your_proc_dir/graynet/spm_cat_gmdensity_cat_aal_range_0_1_nbins25/``



.. _roi_stats:

Computing ROI-wise statistics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When available, ``graynet`` will help compute ROI-wise statistics (individual, not pair-wise) for visualization (e.g. mean and variance of GM density within amygdala), as well as to serve as a baseline for network-level features. Use the ``-r`` or ``--roi_stats`` flag to achieve it. Only one type of processing (ROI stats, or network-level features) can be done at a time.