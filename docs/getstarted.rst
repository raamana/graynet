---------------
Getting Started
---------------

The following steps should help you get started and going quickly:

    - ensure Freesurfer processing is complete and
    - you ran ``recon-all`` with -qcache flag choosing atleast one FWHM value (10 is the default). If not, it is quick to rerun on existing Freesurfer processing.
    - Install graynet using this command in a terminal:  ``pip install graynet``
    - User `script <https://github.com/raamana/graynet/blob/master/scripts/generate_hpc_jobs.py>`_ to produce jobs for processing on a cluster, or even locally.


After processing is done on the cluster and assuming no other errors, you should see a CSV file for each subject, with ``n*(n-1)/2`` values corresponding to the upper-triangular part of the connectivity matrix, where n is the number of ROIs in a atlas you chose (68 for fsaverage or 360 for Glasser2016).

The output folder will be graynet within ``proc_dir`` you choose in the script - determined by this statement in the above script:

.. code-block:: python

    out_dir = pjoin(proc_dir, 'graynet', '{}_{}_fwhm{}_range{}_{}_nbins{}'.format(base_feature, atlas, fwhm, edge_range[0], edge_range[1], num_bins))


which after expansion looks something like based on your choices: ``/your_proc_dir/graynet/freesurfer_thickness_GLASSER2016_fwhm10_range0_5_nbins25/``


If you run into any issues, or have a feature you would like, please let me know here `by opening an issue <https://github.com/raamana/graynet/issues/new>`_.

Thanks for trying out graynet.


Citation
--------

If you found it useful in your research, I'd appreciate if you could cite it using this info:

Pradeep Reddy Raamana. (2017, September 26). *graynet: Subject-wise networks from structural shape and morphological features*. Zenodo. http://doi.org/10.5281/zenodo.997358