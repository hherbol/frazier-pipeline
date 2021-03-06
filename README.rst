Frazier PipeLine
==============================

Pipeline for submitting solubility simulations.  Bonus output may include Unsaturated Mayer Bond Order amongst other things.

Installing
------------------------------

Currently installation involves first installing clancelot:

.. code-block:: bash

   [user@local]~% cd ~; git clone git@github.com:clancylab/clancelot2.0.git

And then cloning this project:

.. code-block:: bash

   [user@local]~% cd ~; git clone git@github.com:hherbol/frazier-pipeline.git

Note, you'll also have to append the frazier-pipeline/pys folder to your PYTHONPATH variable.

.. code-block:: bash

   [user@local]~% echo '\nexport PYTHONPATH="/PATH/TO/FRAZIER/PIPELINE/pys:$PYTHONPATH"' >> ~/.zshrc

Documentation
------------------------------

Documentation is necessary, and the following steps MUST be followed during contribution of new code:

**Setup**

1. Download Sphinx_.  This can be done simply if you have pip_ installed via `pip install -U Sphinx`

2. Wherever you have *frazier-pipeline* installed, you want another folder called *frazier-pipeline-docs* (NOT as a subfolder of frazier-pipeline).

.. code-block:: bash

   [user@local]~% cd ~; mkdir frazier-pipeline-docs; cd frazier-pipeline-docs; git clone -b gh-pages git@github.com:hherbol/frazier-pipeline.git html

3. Forever more just ignore that directory (don't delete it though)

**Adding Documentation**

Documentation is done using ReStructuredText_ format docstrings, the Sphinx_ python package, and indices with autodoc extensions.  To add more documentation, first add the file to be included in `docs/source/conf.py` under `os.path.abspath('example/dir/to/script.py')`.  Secondly, ensure that you have proper docstrings in the python file, and finally run `make full` to re-generate the documentation and commit it to your local branch, as well as the git *gh-pages* branch.

For anymore information on documentation, here is a link to the followed sphinx_tutorial_.

.. _tutorial: https://www.atlassian.com/git/tutorials/using-branches/git-branch
.. _Sphinx: http://www.sphinx-doc.org/en/stable/
.. _pip: https://pip.pypa.io/en/stable/installing/
.. _ReStructuredText: http://docutils.sourceforge.net/docs/user/rst/quickref.html
.. _sphinx_tutorial: https://daler.github.io/sphinxdoc-test/includeme.html

Using FPL
------------------------------

.. automodule:: fpl_auto
    :members:
