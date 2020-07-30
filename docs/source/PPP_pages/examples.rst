==================
Usage and Examples
==================

###############################
Command-line Usage and Examples
###############################

PPP functions may be called at the command-line as shown in this example:

.. code-block:: bash
        
   vcf_filter.py --vcf examples/files/merged_chr1_10000.vcf.gz --filter-only-biallelic --out-format bcf

Details on the usage and arguments of each function may be found within the relevant documentation. In addition, all files specified within these examples may be found within the **examples/files** directory of the PPP repository.

#########################
Module Usage and Examples 
#########################

PPP functions may also be imported from the **pgpipe** module for use within a python script or a `Jupyter Notebook <https://jupyter.org/>`_ as shown in this example:

.. code-block:: python
		
   import pgpipe.vcf_filter as vcf_filter
        
   vcf_filter.run(vcf = 'examples/files/merged_chr1_10000.vcf.gz', filter_only_biallelic = True, out_format = 'bcf')

In comparison to calling functions at the command-line, imported functions require:

* The pgpipe module must be imported
* Each function is called using .run(): **vcf_filter.run() or pgpipe.vcf_filter.run()**
* The use of underscores within arguments rather than dashes: **--out-format vs. out_format**
* Setting the value to True when arguments do not require a value: **--filter-only-biallelic vs. filter_only_biallelic = True**

*************************
Example Jupyter Notebooks
*************************

We have included two example notebooks within the **examples/jupyter** directory of the PPP repository.

.. toctree::
   :maxdepth: 1

   jupyter/example_pipeline_pan.ipynb

.. only:: html

   Jupyter Notebooks may also be download using the following links:

   * :download:`Example Jupyter Pipleine <jupyter/example_pipeline_pan.ipynb>`.
