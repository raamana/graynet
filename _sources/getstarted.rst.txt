---------------
Getting started
---------------

``graynet`` currently offers the following streams of processing:


  - :doc:`cortical` : using vertex-wise ROIs defined on the cortex. This is useful to analyze network-level features based on cortical thickness, curvature, sulcal depth and gyrification. Base-level features can easily obtained from running `Freesurfer <https://surfer.nmr.mgh.harvard.edu/>`_
  - :doc:`volumetric` stream : using voxel-wise ROIs defined over the whole brain relying on a volumetric atlas. This is useful to analyze network-level features based on gray matter density dervied from voxel-based morphometry (VBM) or similar approaches. Base-level features can easily obtained from `CAT12 toolbox <http://www.neuro.uni-jena.de/cat/>`_ within the SPM ecosystem.

We plan to offer the following soon:


  - support for other VBM tools such as `SPM <http://www.fil.ion.ucl.ac.uk/spm/software/spm12/>`_, `ANTs <http://stnava.github.io/ANTs/>`_, `FSL <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLVBM>`_, `Freesurfer <https://surfer.nmr.mgh.harvard.edu/>`_, or another suitable package, within the :doc:`volumetric` stream, .
  - Support for additional input formats for the :doc:`cortical` stream: `CIVET <http://www.bic.mni.mcgill.ca/ServicesSoftware/BasicUsageOfCIVET>`_ , `ANTs <http://stnava.github.io/ANTs/>`_ etc.
  - If you are interested in contributing, please take a look at the :doc:`contributing` document and reach out to me. Thanks!

In both streams, in addition to the computation of pair-wise network-level features, ``graynet`` will help you compute ROI-wise statistics (individual, not pair-wise) for visualization (median thickness in PCG, or variance in GM density within amygdala), as well as to serve as a baseline for network-level features.


The package offers both :doc:`usage_cli` and an :doc:`API`, to better integrate with your workflow. However, the CLI is the recommended/most-tested gateway.

For uniform processing across subjects, ``graynet`` needs:

  - an atlas with pre-defined ROIs
  - each subject be registered to atlas (vertex- or voxel-wise), so ROIs correspond across all subjects.
  - the extracted features are in a format readable by ``nibabel``

The following steps should help you get started and going quickly:

    - Install graynet using this command in a terminal:  ``pip install -U graynet``
    - Refer to :doc:`cortical` and :doc:`volumetric` pages for more details on the individual streams of processing and examples.


If you run into any issues, or have a feature wish or suggestions, please let me know here `by opening an issue <https://github.com/raamana/graynet/issues/new>`_.

Thanks for trying out graynet. **I'd appreciate if you can cite it** using the details in :doc:`citation`.

