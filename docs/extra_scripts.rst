
Convenience/useful scripts
--------------------------------

To automate the the process of generating the above job files and their submission for processing on HPC clusters, I have included a `script <https://github.com/raamana/graynet/blob/master/scripts/generate_hpc_jobs.py>`_ . Please read the instructions inside and run it.

There are few other convenience scripts in the `scripts <https://github.com/raamana/graynet/blob/master/scripts/>`_ folder - take a look and feel free to modify for your tasks.


We encourage you to adopt `GraphML` format to store and analyze the networks extrated from graynet, however, there maybe simple use cases you prefer simpler CSVs. To convert `GraphML` files to CSV format (contianing just weight values, and nothing else), use the script https://github.com/raamana/graynet/blob/master/scripts/convert_graphml_to_csv.py
