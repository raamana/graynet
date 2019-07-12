------------
Installation
------------

This package could easily be installed via:

.. code-block::

    pip install -U graynet


``graynet`` is written in Python completely and is distributed via the Python Package Index. It requires python (Free and awesome!) to be installed and `pip` to install. You most likley have them already.


If you do not have python or pip installed, please follow the links below for instructions on installing them:

 - Python 3 or higher: https://www.python.org/downloads/ (pip comes packaged with Python)
 - **Windows** users:
 
   - Follow this guide https://www.ics.uci.edu/~pattis/common/handouts/pythoneclipsejava/python.html or https://www.python.org/downloads/windows/
   - and read these FAQ to familiarize themselves with typical questions using these FAQ: https://docs.python.org/3/faq/windows.html


.. hint::

    Note as this package is mainly geared towards batch processing (typically done on high-performance computing clusters). Hence, this has been tested only on Linux, and **has not been tested on Windows**. The API should work similarly across all OSes as python is cross-platform by design. However, the command line interface on Windows might be more difficult to use compared to Linux (because of the horrible nature of path specification and handling on Windows, and lack of powerful Terminals for easy file management).


Requirements
------------

 - numpy
 - nibabel
 - hiwenet


Citation
--------

If you found any parts of graynet to be useful in your research, I'd appreciate if you could cite the following papers noted at :doc:`citation`.