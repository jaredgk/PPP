========
Examples
========

#################
Function Examples
#################

PPP functions may be called at the command-line as shown in this example:

.. code-block:: bash
        
        vcf_filter.py --vcf examples/files/merged_chr1_10000.vcf.gz --filter-only-biallelic --out-format bcf

Details on the usage of each specific function may be found within the *Example usage* section of the function’s documentation. In addition, all files shown within these examples may be found within **examples/files** directory of the PPP repository.

##########################
Jupyter Notebook Pipelines
##########################

All PPP functions may also be used within a `Jupyter Notebook <https://jupyter.org/>`_. We have included two example notebooks.

.. toctree::
   :maxdepth: 1

   jupyter/example_pipeline_pan.ipynb

.. only:: html

	The Jupyter Notebooks may also be download:

   * :download:`Example Jupyter Pipleine <jupyter/example_pipeline_pan.ipynb>`.