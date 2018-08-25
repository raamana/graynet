------------
Installation
------------

This package could easily be installed via:

``pip install -U graynet``

``graynet`` is written in Python completely and is distributed via the Python Package Index. It requires python (Free and awsome!) to be installed and `pip` to install. You most likley have them already. If you do not have python or pip installed, please follow the links below for instructions on installing them:

 - Python 3 or higher: https://www.python.org/downloads/ (pip comes packaged with Python)
 - Windows users:
 
   - can follow this guide https://www.ics.uci.edu/~pattis/common/handouts/pythoneclipsejava/python.html or https://www.python.org/downloads/windows/
   - and read these FAQ to familiarize themselves with typical questions using these FAQ: https://docs.python.org/3/faq/windows.html
 - Note as this package is mainly geared towards batch processing (typically done on high-performance computing clusters). Hence, this has been tested only on Linux, and has not been tested on Windows. The API should work similarly, but the command line interface on Windows might be slightly less easy compared to Windows.


Requirements
------------

 - numpy
 - nibabel
 - hiwenet


Citation
--------

If you found any parts of graynet to be useful in your research, I'd appreciate if you could cite one or more of the following papers in that order:

 - Raamana, P. R. & Strother, S. C., (2018). *graynet: single-subject morphometric networks for neuroscience connectivity applications* (Version 0.3.9). Zenodo. http://doi.org/10.5281/zenodo.1403304
 - Pradeep Reddy Raamana. *graynet: Subject-wise networks from structural shape and morphological features*. Zenodo. http://doi.org/10.5281/zenodo.997358
 - Raamana, P.R. and Strother, S.C., 2017, Impact of spatial scale and edge weight on predictive power of cortical thickness networks bioRxiv 170381 http://www.biorxiv.org/content/early/2017/07/31/170381.
 - Raamana, P. R., Weiner, M. W., Wang, L., Beg, M. F., & Alzheimer's Disease Neuroimaging Initiative. (2015). Thickness network features for prognostic applications in dementia. Neurobiology of aging, 36, S91-S102.