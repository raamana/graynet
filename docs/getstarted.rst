---------------
Getting started
---------------

The following steps should help you get started and going quickly:

    - Ensure Freesurfer processing is complete and
    - Ensure you ran ``recon-all`` with ``-qcache`` flag choosing atleast one FWHM value (10 is the default), for the feature you are interested in analyzing (e.g. thickness). If not already run, it is quick to rerun on existing Freesurfer processing.
    - Install graynet using this command in a terminal:  ``pip install -U graynet``
    - Use this `script <https://github.com/raamana/graynet/blob/master/scripts/generate_hpc_jobs.py>`_ to produce jobs for processing on a cluster, or even locally.


After processing is done on the cluster and assuming no other errors, you should see a CSV file for each subject, with ``n*(n-1)/2`` values corresponding to the upper-triangular part of the connectivity matrix, where n is the number of ROIs in a atlas you chose (68 for fsaverage or 360 for Glasser2016).

The output folder will be graynet within ``proc_dir`` you choose in the script - determined by this statement in the above script:

.. code-block:: python

    out_dir = pjoin(proc_dir, 'graynet', '{}_{}_fwhm{}_range{}_{}_nbins{}'.format(base_feature, atlas, fwhm, edge_range[0], edge_range[1], num_bins))


which after expansion looks something like based on your choices: ``/your_proc_dir/graynet/freesurfer_thickness_GLASSER2016_fwhm10_range0_5_nbins25/``


Qcache recon-all flag
----------------------------------

The Qcache recon-all flag does the following:

    - resample data (e.g. thickness, curv, sulc) onto the average subject (called fsaverage)
    - smooth it at a range of FWHM (full-width/half-max) values, usually 0, 5, 10, 15, 20, and 25mm.

Check https://surfer.nmr.mgh.harvard.edu/fswiki/qcache for more info

The parcellation for the `Glasser et al 2016 <http://www.nature.com/nature/journal/vaop/ncurrent/full/nature18933.html>`_ atlas has been resampled onto the fsaverage space already, courtesy of `Kathryn Mills <https://figshare.com/articles/HCP-MMP1_0_projected_on_fsaverage/3498446>`_.


If you run into any issues, or have a feature you would like, please let me know here `by opening an issue <https://github.com/raamana/graynet/issues/new>`_.

Thanks for trying out graynet.


Citation
--------

If you found it useful in your research, I'd appreciate if you could cite it using this info:

Pradeep Reddy Raamana. (2017, September 26). *graynet: Subject-wise networks from structural shape and morphological features*. Zenodo. http://doi.org/10.5281/zenodo.997358