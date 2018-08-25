---
title: 'single-subject morphometric networks for neuroscience connectivity applications'
tags:
  - neuroscience
  - network
  - morphometry
  - connectivity
  - gray matter
  - graph
  - histogram
  - freesurfer
authors:
 - name: Pradeep Reddy Raamana
   orcid: 0000-0003-4662-0558
   affiliation: 1
 - name: Stephen C. Strother
   orcid: 0000-0002-3198-217X
   affiliation: 1, 2
affiliations:
 - name: Rotman Research Institute, Baycrest Health Sciences, Toronto, ON, Canada
   index: 1
 - name: Department of Medical Biophysics, University of Toronto, Toronto, ON, Canada
   index: 2
date: 23 August 2018
doi: 10.5281/zenodo.1403304
bibliography: paper.bib
---

# Summary

Connectivity, and network-level features in general, have proven to be valuable tools in several aspects of neuroscience research. Although network analysis is rooted in analysis of functional MRI data, it has recently gained traction in the analyses of morphometric features such as cortical thickness [@evans2013networks]. Such networks of anatomical covariance (derived based on distributions of features across a group of subjects) provided insight into changes caused by various brain disorders. When we individualize this approach to enable extraction of single-subject network features, they provided several useful insights [@tijms2012similarity,@raamana2015thickness,@palaniyappan2015abnormalities,@xu2017abnormalities]. Moreover, network-level features demonstrated potential for prognostic applications [@raamana2015thickness,@raamana2014novel], in addition to being robust to changes in scale and edge weight metrics [@Raamana2017biorxiv].

However, deriving these network-level features from input T1w-MRI data is non-trivial. With this fully-open-source and pure-python library `graynet`, we attempt to make this task relatively much easier and extend it support all currently available morphometric features. Currently, it interfaces directly with the outputs produced by Freesurfer [@fischl2012freesurfer] and supports vertex-wise data. We plan to extend this to support volumetric atlases and voxel-wise features such as gray matter density. Together with many convenience scripts (e.g. to launch jobs on the high-performance cluster and assemble the outputs produced), we believe `graynet` makes an useful addition to the neuroimaging in python open source ecosystem.

`graynet` is dependent on the following libraries: `nibabel` [@brett2016nibabel], `networkx` [@hunter2007matplotlib], `numpy` [@oliphant2007python,@walt2011numpy] and `hiwenet` [@raamana2017hiwenet].

# Acknowledgement

Pradeep Reddy Raamana is grateful for the support of the Canadian Biomarker Integration Network for Depression (CAN-BIND) and Ontario Neurodegenerative Disease Research Initiative (ONDRI), which are two integrated discovery programs of the Ontario Brain Institute (OBI), Canada. OBI is an independent non-profit corporation, funded partially by the Ontario government. The opinions, results, and conclusions are those of the authors and no endorsement by the OBI is intended or should be inferred.

# References
