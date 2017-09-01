# graynet

Individualized single-subject (covariance) networks from T1 mri features such as cortical thickness, gray matter density, subcortical morphometric features, gyrification and curvature. 

Applicable for biomarker development for the detection of many brain disorders such as Alzheimer's, Parkinson's, Schizophrenia and the like - see below.

## References
The following publications outline typical use cases:
 * Raamana, P.R. and Strother, S.C., 2017, Impact of spatial scale and edge weight on predictive power of cortical thickness networks bioRxiv 170381 http://www.biorxiv.org/content/early/2017/07/31/170381. doi: https://doi.org/10.1101/170381
 * Tijms, B. M., Seriès, P., Willshaw, D. J., & Lawrie, S. M. (2012). Similarity-based extraction of individual networks from gray matter MRI scans. Cerebral Cortex, 22(7), 1530-1541.
 * Palaniyappan, L., Park, B., Balain, V., Dangi, R., & Liddle, P. (2015). Abnormalities in structural covariance of cortical gyrification in schizophrenia. Brain Structure and Function, 220(4), 2059-2071.
 * Xu, J., Zhang, J., Zhang, J., Wang, Y., Zhang, Y., Wang, J., ... & Zhang, Y. (2017). Abnormalities in structural covariance of cortical gyrification in Parkinson's disease. Frontiers in Neuroanatomy, 11.
 * Raamana, P. R., Wen, W., Kochan, N. A., Brodaty, H., Sachdev, P. S., Wang, L., & Beg, M. F. (2014). Novel ThickNet features for the discrimination of amnestic MCI subtypes. NeuroImage: Clinical, 6, 284-295.

## Installation

`pip install graynet`

## Usage

A rough example of usage can be:

```python
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

```
