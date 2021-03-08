============
Installation
============

#########
From PyPi
#########

The PPP can also be easily installed via the PyPi repository via pip:

.. code-block:: bash

    pip install py-popgen
    
##########   
From Conda
##########

The PPP has conda packages available for python versions 3.6 and 3.7. To install in a clean environment, run the following:

.. code-block:: bash

    conda create -n py-popgen python=3.7.7
    conda activate py-popgen
    conda install -c jaredgk -c bioconda py-popgen

###########
From Source
###########

The most current version of the PPP can be installed by obtaining the source code from the PPP GitHub repository. This can be done with:

.. code-block:: bash

    git clone https://github.com/jaredgk/PPP

To install the local repository copy and allow edits to the source code to be included with imports without any additional steps, run the following commands:

.. code-block:: bash
    
    cd PPP
    pip install -e . 
    
To install the repository without pip, run the following (note that any modifications to the source code will not be used at runtime unless the setup command is run again):

.. code-block:: bash

    cd PPP
    python setup.py install 

############
Dependencies
############

If installing PPP from source, multiple python and non-python depencencies must also be installed. 

-------------------
Python Dependencies
-------------------

The PPP requries a number of python libraries, including:

* `The SciPy Ecosystem <https://www.scipy.org/about.html>`_ (i.e. numpy, scipy, pandas, matplotlib, etc.)
* `Pysam <https://github.com/pysam-developers/pysam>`_
* `Biopython <https://biopython.org/>`_  
* `Cython <https://cython.org/>`_  
* `rpy2 <https://rpy2.readthedocs.io/>`_

We recommend users install and maintain these libraries using either `pip <https://pypi.org/project/pip/>`_ or `Anaconda 3 <https://www.anaconda.com/distribution/#download-section>`_.

------------------
Other Dependencies
------------------

The PPP also requries a number of executables to be installed, including:

* `VCFtools <https://vcftools.github.io/index.html>`_
* `BCFtools, Samtools, and HTSlib <http://www.htslib.org/>`_
* `plink 1.9 <https://www.cog-genomics.org/plink2/>`_
* `plink 2.0 <https://www.cog-genomics.org/plink/2.0/>`_
* `SHAPEIT <https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html>`_
* `Beagle 5.0 <https://faculty.washington.edu/browning/beagle/beagle.html>`_
* `Picard <https://broadinstitute.github.io/picard/>`_

Please note that VCFtools, BCFtools, Samtools, HTSlib, plink 1.9, plink 2.0, and SHAPEIT may be installed using `Anaconda 3 <https://www.anaconda.com/distribution/#download-section>`_.


