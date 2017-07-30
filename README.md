# graynet
Individualized single-subject networks from T1 mri features such as cortical thickness and gray matter density. 


## References
A publication outlining one use case is here:
[Raamana, P.R. and Strother, S.C., 2016, June. Novel histogram-weighted cortical thickness networks and a multi-scale analysis of predictive power in Alzheimer's disease. In Pattern Recognition in Neuroimaging (PRNI), 2016 International Workshop on (pp. 1-4). IEEE.](http://ieeexplore.ieee.org/abstract/document/7552334/)

Another poster describing it can be found here: https://doi.org/10.6084/m9.figshare.5241616

## Installation

`pip install -U graynet`

## Usage

A rough example of usage can be:

```python
from graynet import extract as graynet
import os

base_project_dir = '/project/mydataset'
freesurfer_dir = os.path.join(base_project_dir, 'freesurfer_v5.3')
base_feature = 'thickness'
fwhm = 10
parcellation = 'nearestneighbour'
patch_size = 1000 # number of vertices per patch in this case

out_dir = os.path.join(base_project_dir, 'graynet', base_feature)

graynet(freesurfer_dir, base_feature, parcellation, fwhm=fwhm, size = patch_size)

```
