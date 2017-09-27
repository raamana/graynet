------------------
Examples using API
------------------

A rough example of usage when using the graynet API is shown below, although I strongly recommend using the command line interface given it is geared towards batch processing on HPC clusters and there is a `script <https://github.com/raamana/graynet/blob/master/scripts/generate_hpc_jobs.py>`_ to help you with large amount of processing.

.. code-block:: python

    import graynet
    import os

    base_project_dir = '/project/mydataset'
    input_dir = os.path.join(base_project_dir, 'freesurfer_v5.3')
    base_feature = 'thickness' # 'gmdensity', 'curv' or 'sulc' or 'gyrification'
    fwhm = 10
    atlas = 'GLASSER2016'
    patch_size = 'default' # number of vertices per patch in this case

    out_dir = os.path.join(base_project_dir, 'graynet', base_feature)

    graynet.extract(input_dir, base_feature, atlas, smoothing_param = fwhm, size = patch_size)



**TODO** add more examples