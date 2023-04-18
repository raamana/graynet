# graynet

[![DOI](http://joss.theoj.org/papers/10.21105/joss.00924/status.svg)](https://doi.org/10.21105/joss.00924)
[![saythanks](https://img.shields.io/badge/say-thanks-ff69b4.svg)](https://saythanks.io/to/raamana)

# News

 - Volumetric version of graynet is available now - [check here](https://raamana.github.io/graynet/volumetric.html)


# Overview

Individualized single-subject networks from T1-weighted magnetic resonance imaging (MRI) features such as:
  - vertex-wise features such as cortical thickness, gyrification and curavature
  - volumetric features such as gray matter density (T1w images) or metabolic uptake (PET) or another voxel-wise feature
  - Subcortical morphometric features
  - or any other similar feature distributed over a domain (cortex, or whole brain) to enable compuatation of inter-regional edge weights

Applicable for whenever network-level features are useful, among which common use cases are: 
 - Biomarker development.
 - Brain-behaviour relationships (e.g. for the diagnosis and prognosis of many brain disorders such as Alzheimer's, Parkinson's, Schizophrenia and the like).
 - Aging (changes in network properties over age and their relations to other variables).

Docs: https://raamana.github.io/graynet/

Quick illustration:

![graynet_flyer](docs/vis/graynet_flyer.jpg)

## Installation

`pip install -U graynet`

## Citation

If you found any parts of graynet to be useful in your research, I'd appreciate if you could cite the software paper in JOSS below, as well as the methods paper that motivated the tool development in that order:

 - Raamana et al., (2018). graynet: single-subject morphometric networks for neuroscience connectivity applications. Journal of Open Source Software, 3(30), 924, https://doi.org/10.21105/joss.00924
 - Raamana, P. R., & Strother, S. C. (2020), “Does size matter? Relationship between predictive power of single subject morphometric networks to spatial scale and edge weight”, Brain Structure and Function, 225(8), 2475-2493. DOI: 10.1007/s00429-020-02136-0
 - Raamana, P. R., Weiner, M. W., Wang, L., Beg, M. F., & Alzheimer's Disease Neuroimaging Initiative. (2015). Thickness network features for prognostic applications in dementia. Neurobiology of aging, 36, S91-S102. https://www.sciencedirect.com/science/article/pii/S0197458014005521

---

[![saythanks](https://img.shields.io/badge/say-thanks-ff69b4.svg)](https://saythanks.io/to/raamana)

