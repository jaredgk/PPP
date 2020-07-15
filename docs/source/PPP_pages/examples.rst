========
Examples
========

#################
Function Examples
#################

PPP functions may be called at the command-line as shown in this example:

.. code-block:: bash
        
        vcf_filter.py --vcf examples/files/merged_chr1_10000.vcf.gz --filter-only-biallelic --out-format bcf

Details on the usage of a specific function may be found within the *Example usage* section of the function in question. In addition, all example files used may be found within **examples/files** directory.  


##########################
Jupyter Notebook Pipelines
##########################

All PPP functions may also be used within a `Jupyter Notebook <https://jupyter.org/>`_. We have included some examples below:

.. toctree::
   :maxdepth: 1
   
   example_pipeline_pan.ipynb