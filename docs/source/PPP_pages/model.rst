=======================
Model File and Creation
=======================

A core aspect of the PPP is the use of Model files, JSON-based files used to assign and store **population models**. A population model primarily consists of: the populations within the model; the individuals in each population; and a population tree. Model files offer various benefits within the PPP: i) automatic assignment of relevant populations, individuals, or other potential meta-data; ii) simplified process to examine multiple models; and iii) a single repository of all relevant meta-data.


Model files may be created and edited using our model creator. 

.. toctree::
   :maxdepth: 1

   Functions/model_creator


An example Model file may be seen below:

.. code-block:: bash

    [
    {
        "name": "2Pop",
        "tree": "(Troglodytes,Verus);",
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
