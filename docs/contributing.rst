Contributing
-------------

Thank you for the checking out ``graynet``, and for your interest in contributing.


Please check the following guidelines for contributing: https://opensource.guide/how-to-contribute/


Developer Guidelines
---------------------

The preferred way to contribute to ``graynet`` is to fork the `main repository <https://github.com/raamana/graynet/>`__ on GitHub,
then submit a "pull request" (PR):

 1. `Create an account <https://github.com/join>`_ on
    GitHub if you do not already have one.

 2. Fork the `project repository <https://github.com/raamana/graynet>`__: click on the 'Fork'
    button near the top of the page. This creates a copy of the code under your
    account on the GitHub server. For more details on how to fork a
    repository see `this guide <https://help.github.com/articles/fork-a-repo/>`_.

 3. Clone this copy to your local disk::

        $ git clone git@github.com:YourLogin/graynet.git

 4. Create a branch to hold your changes::

        $ git checkout -b my-feature

    and start making changes. Never work in the ``master`` branch!

 5. Work on this copy, on your computer, using Git to do the version
    control. When you're done editing, do::

        $ git add modified_files
        $ git commit

    to record your changes in Git, then push them to GitHub with::

        $ git push -u origin my-feature

Finally, follow `these <https://help.github.com/articles/creating-a-pull-request-from-a-fork>`_ instructions to create a pull request from your fork.

.. note::

  In the above setup, your ``origin`` remote repository points to the repository under your github account ``YourLoginName/graynet.git`` .
  If you wish to fetch/merge from the main repository instead of your forked one, you will need to add another remote
  to use instead of ``origin``. If we choose the name ``upstream`` for it, the command will be::

        $ git remote add upstream https://github.com/raamana/graynet.git

If any of the above seems like magic to you, then look up the `Git documentation
<https://git-scm.com/documentation>`_ and the `Git development workflow
<http://docs.scipy.org/doc/numpy/dev/gitwash/development_workflow.html>`_ on the
web.

If some conflicts arise between your branch and the ``master`` branch, you need
to merge ``master``. The command will be::

  $ git merge master

with ``master`` being synchronized with the ``upstream``.

Subsequently, you need to solve the conflicts. You can refer to the `Git
documentation related to resolving merge conflict using the command line
<https://help.github.com/articles/resolving-a-merge-conflict-using-the-command-line/>`_.

.. note::

    These guidelines are adapted from ``scikit-learn`` docs, credit to the original contributors.


Testing
--------


.. code-block:: python

    cd <<path to root folder of graynet repo cloned>>
    cd graynet/tests
    pytest


that will run all tests in `test_graynet.py` and it should report no errors (like ``12 passed`` below):

.. code-block:: bash

    $ 08:38:39 SQuark-3 tests >>  pwd
    /Users/Reddy/dev/graynet/graynet/tests
    $ 08:38:48 SQuark-3 tests >>  pytest
    ===================================================================================== test session starts ======================================================================================
    platform darwin -- Python 3.6.2, pytest-3.4.2, py-1.5.2, pluggy-0.6.0
    rootdir: /Users/Reddy/dev/graynet, inifile:
    plugins: hypothesis-3.30.1
    collected 12 items

    test_graynet.py ............                                                                                                                                                             [100%]

    ======================================================================================= warnings summary =======================================================================================
    graynet/tests/test_graynet.py::test_run_roi_stats_via_API
      /Users/Reddy/anaconda/envs/py36/lib/python3.6/site-packages/scipy/stats/stats.py:2831: FutureWarning:

      Using a non-tuple sequence for multidimensional indexing is deprecated; use `arr[tuple(seq)]` instead of `arr[seq]`. In the future this will be interpreted as an array index, `arr[np.array(seq)]`, which will result either in an error or a different result.


    -- Docs: http://doc.pytest.org/en/latest/warnings.html
    =========================================================================== 12 passed, 1 warnings in 106.92 seconds ============================================================================


If you would like to add more tests, add a new function in ``test_graynet.py``, or add a new file to ``graynet/tests`` folder with a filename starting with ``test_``.

For more info, check the following links

 - on testing python packages : https://www.python-course.eu/python3_tests.php
 - ``pytest`` framework : https://docs.pytest.org/en/latest/contents.html
 - pytest tutorial : https://semaphoreci.com/community/tutorials/testing-python-applications-with-pytest
 - Good Integration Practices from pytest devs : https://docs.pytest.org/en/latest/goodpractices.html
 - comprehensive coverage: https://wiki.python.org/moin/PythonTestingToolsTaxonomy