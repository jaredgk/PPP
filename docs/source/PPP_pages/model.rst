==========
Model File
==========

A unique benefit of the PPP is the ability to use the Model files. Model files are JSON-based files used to store population models and their relevant details, such as: the populations within the model; the individuals in each population; the population tree; and other potential meta-data as needed. Model files allow the PPP functions to automatically assign various parameters and serves as the repository for all model-related information.

.. toctree::
   :maxdepth: 1

   Functions/model_creator


An example Model file may be seen below:

.. code-block:: bash

    [
    {
        "name": "2Pop",
        "pops": {
            "Verus": {
                "inds": [
                    "Pan_troglodytes_verus-9668_Bosco",
		    "Pan_troglodytes_verus-9730_Donald",
		    "Pan_troglodytes_verus-A956_Jimmie",
		    "Pan_troglodytes_verus-Clint",
		    "Pan_troglodytes_verus-X00100_Koby"
                ]
            },
            "Troglodytes": {
                "inds": [
                    "Pan_troglodytes_troglodytes-A957_Vaillant",
                    "Pan_troglodytes_troglodytes-A958_Doris",
                    "Pan_troglodytes_troglodytes-A959_Julie",
                    "Pan_troglodytes_troglodytes-A960_Clara"
                ]
            }
        }
    }
    ]
