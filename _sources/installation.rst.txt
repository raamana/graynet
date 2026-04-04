------------
Installation
------------

Install graynet from PyPI:

``pip install -U graynet``

graynet 2.0 requires Python 3.10 or newer.

Core dependencies are installed automatically. The main runtime dependencies are
NumPy, NiBabel, NetworkX, PyArrow, and hiwenet.

Quick sanity check
------------------

After installation, the CLI should be available:

.. code-block:: bash

    graynet --version

You can also verify the top-level API:

.. code-block:: python

    from graynet import extract, extract_multiedge, roiwise_stats_indiv
    from graynet import load_run, get_edge_values, export_to_nx


Requirements
------------

- Python 3.10+
- NumPy
- NiBabel
- PyArrow
- hiwenet


Citation
--------

If you found any parts of graynet to be useful in your research, I'd appreciate if you could cite the following papers noted at :doc:`citation`.
