Cortical graynet
-------------------

For the cortical stream of processing, ``graynet`` relies on outputs from Freesurfer processing. If you processed your data (to produce thickness values), just resample them on fsaverage atlas and store it in a format readable by ``nibabel`` package.

The following steps should help you get started and going quickly:

    - Ensure Freesurfer processing is complete. Check :ref:`run_fs` for more details
    - Ensure you ran ``recon-all`` with ``-qcache`` flag choosing atleast one FWHM value (10 is the default). If not already run, it is quick to rerun on existing Freesurfer processing. See :ref:`qcache_flag` for more details
    - It is **highly** recommended that you perform quality control on the Freesurfer outputs before you use them for analysis. Such QC is quite easy now using `visualqc <https://raamana.github.io/visualqc/readme.html>`_.

If graynet is installed without error, you should have ``graynet`` command in your path. Just type ``graynet`` and you should see it display options and usage instructions, as presented in :doc:`usage_cli` page.

Suppose say

  - you already have Freesurfer processing done at ``/work/project/freesurfer_reconall``
  - for 5 subjects with identifiers are stored in subject_ids.txt (one per line)
  - you want to analyze network features based on cortical thickness values, which have been resampled on fsaverage surface (which is the default) at smoothing level of FWHM 10 (default)
  - use the ``manhattan`` method to compute histogram distance between two cortical ROIs (assuming all subjects resampled to same atlas)
  - and store the networks extracted in /work/project/graynet

You could achieve it all with a single command:

.. code-block:: bash

    cd /work/project
    graynet -s subject_ids.txt -f freesurfer_thickness -i /work/project/freesurfer_reconall \
        -w manhattan -a fsaverage -p 10 -o /work/project/graynet

That's it! By the time you can get your coffee or stretch your legs, you should have graynet processing done.

Suppose, you prefer to analyze ROIs as defined by a multimodal parcellation published in Glasser 2016 atlas and also curious to try multiple weight metrics -  ``chebyshev`` and ``cosine``, you can achieve it via the following command:

.. code-block:: bash

    cd /work/project
    graynet -s subject_ids.txt -f freesurfer_thickness -i /work/project/freesurfer_reconall \
        -w manhattan chebyshev cosine -a Glasser2016 -p 10 -o /work/project/graynet


You could also study curvature and sulcal depth features by simply adding more features to the ``-f``, such as ``freesurfer_curv`` and ``freesurfer_sulc``.


Using a different Atlas
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

when you want to use a different atlas, you must have already processed and resampled the desired features in that atlas space, with a separate ``mgh`` file for each atlas. For example, when you run Freesurfer's ``QCache`` with ``fsaverage`` atlas, you would see the files with the following pattern  ``?h.thickness.fwhm10.fsaverage.mgh`` in the ``surf`` folder of the subject's Freesurfer output (e.g. in ``$SUBJECTS_DIR/$SUBID``).

If you want to use a different atlas such as ``glasser2016``, graynet expects to see files of the pattern ``?h.thickness.fwhm10.glasser2016.mgh``. The ``glasser2016`` atlas distributed with graynet is already in vertex-wise correspondence with ``fsaverage``, so you can simply add symbolic links (with ``symlink`` command in Unix variants of the OS) for ``?h.thickness.fwhm10.glasser2016.mgh`` to point to ``?h.thickness.fwhm10.fsaverage.mgh`` as a workaround.

That said, if you use a slightly different version of Glasser atlas or something else say ``XYZ42`` e.g. by pointing ``--atlas`` argument to a Freesurfer folder (say ``/work/atlases/XYZ42/freesurfer_parc``) that may not be in vertex-wise correspondence with fsaverage, you need to run the Qcache properly and produce the corresponding features files resampled on the space of ``XYZ42`` and the files must be named accordingly: ``?h.thickness.fwhm10.XYZ42.mgh``. Note I used fwhm as 10 here for illustrative purposes, it could be anything else you choose.


Processing large datasets or atlases with many ROIs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

However, if you have 100s of subjects or your atlas has a large number of ROIs (like the Glasser 2016 which has 360), this computation can take a long time. It is advisable to move the processing to a larger server or computer cluster, by delegating few subjects to a single job. You could leverage this `script <https://github.com/raamana/graynet/blob/master/scripts/generate_hpc_jobs.py>`_ to process any large dataset at once (produce jobs for processing on a cluster).

    - If you have access to a computer cluster, you can submit the generated jobs to the cluster.
    - If you don't have access to a cluster, you could still use the generated job/script files to run them locally - few at a time.

After processing is done (on the cluster or locally), and assuming no other errors, you should see a CSV file for each subject, with ``n*(n-1)/2`` values corresponding to the upper-triangular part of the connectivity matrix, where n is the number of ROIs in a atlas you chose (68 for fsaverage or 360 for Glasser2016).

The output folder will be graynet within ``proc_dir`` you choose in the script - determined by this statement in the above script:

.. code-block:: python

    out_dir = pjoin(proc_dir, 'graynet', '{}_{}_fwhm{}_range{}_{}_nbins{}'.format(base_feature, atlas, fwhm,
                                                edge_range[0], edge_range[1], num_bins))


which after expansion looks something like based on your choices: ``/your_proc_dir/graynet/freesurfer_thickness_GLASSER2016_fwhm10_range0_5_nbins25/``


Take a look at :doc:`extra_scripts` also.

.. _roi_stats_ctx:

Computing ROI-wise statistics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``graynet`` also helps you compute ROI-wise statistics (individual, not pair-wise) for visualization (median thickness in PCG), as well as to serve as a baseline for network-level features. Use the ``-r`` or ``--roi_stats`` flag to achieve it. Only one type of processing (ROI stats, or network-level features) can be done at a time.


.. _run_fs:

How to run Freesurfer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you are new to Freesurfer, please:

    - follow this `beginners guide <https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferBeginnersGuide>`_
    - leverage various `options <https://surfer.nmr.mgh.harvard.edu/fswiki/recon-all>`_ available
    - include ``-qcache`` flag choosing atleast one FWHM value (10 is the default) - see :ref:`qcache_flag`.
    - process the subjects till "recon-all is finished without error"

In typical scenarios (when T1 mri scans do not *special* processing to handle any artefacts), running Freesurfer would boil down to running command:

.. code-block:: bash

    recon-all -all -sd /project/processed/freesurfer -subjid ABC_0001 -i /project/raw/ABC_0001/mri.nii -qcache

.. _qcache_flag:

Qcache recon-all flag
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Qcache recon-all flag does the following:

    - resample data (e.g. thickness, curv, sulc) onto the average subject (called fsaverage)
    - smooth it at a range of FWHM (full-width/half-max) values, usually 0, 5, 10, 15, 20, and 25mm.

We encourage the use of default behaviour (simple inclusion of ``-qcache``), which will smooth the surface data for thickness, curv, sulc, area and jacobian_white at 0, 5, 10, 15, 20, and 25 mm FWHM.

If you already have a completed run of Freesurfer, you can use the following command to run the ``-qcache``

.. code-block:: bash

    recon-all -s $SUBJECT_ID -qcache

which will run it with default parameters i.e. all measures available such as thickness, curvature, sulcal depth, area and the like. In addition, the default behaviour is to smooth them at 5mm FWHM to 25mm FWHM in steps of 5mm. You can control for a fewer combinations (e.g. one measure and one smoothing level) with:

.. code-block:: bash

    recon-all -s $SUBJECT_ID -qcache -measure thickness -fwhm 10


If you are interested in resampling the data on to a different atlas or process features outside Freesurfer structure, check https://surfer.nmr.mgh.harvard.edu/fswiki/qcache for more info.

The parcellation for the `Glasser et al 2016 <http://www.nature.com/nature/journal/vaop/ncurrent/full/nature18933.html>`_ atlas has been resampled onto the fsaverage space already, courtesy of `Kathryn Mills <https://figshare.com/articles/HCP-MMP1_0_projected_on_fsaverage/3498446>`_.

If you run into any issues, or have a feature you would like, please let me know here `by opening an issue <https://github.com/raamana/graynet/issues/new>`_.
