===========
Development
===========

As it would be unlikely that the PPP would encompass all desired procedures/methods, we has included this section to aid in the development of additional functionality for the PPP.

######################
Development Guidelines 
######################

**In general:**

* Functions developed for the PPP should be modular – i.e. function independently – primarily to maintain a flexible platform, this is especially important when considering `Galaxy <https://galaxyproject.org/>`_. integration.
* Functions should be able to be imported from **pgpipe**, so that multiple functions may be used within a single python script or `Jupyter notebook <https://jupyter.org/>`_.
* Functions should support the use of Model files for automatic assignment of relevant meta-data. *Note: details on the model file class may be found below*
* Functions that use third-party software should include the relevant reference(s) within log files. 

#################
Using PPP Classes
#################

*********
VCF Class
*********

*italic*
**BOLD**
outside link: `Google <https://google.com//>`_.
internal link: :ref:`examples`

.. code-block:: bash
        
   my_code.py

***************
ModelFile Class
***************
Model files may be read using the function **read_model_file()** found within the model.py module. The reader accepts the filename (as a string) of the model file and returns a **ModelFile** class object as shown below:

.. code-block:: python

   # Read in the model file
   model_file = pgpipe.model.read_model_file(“path/to/my_model_file.model”)

A **ModelFile** class object primarily behaves as a dictionary of **Model** class objects; the keys of the dictionary are the model names (as strings) of the **Model** class objects whereas the values are the objects themselves. Therefore, a **Model** named *2Pop* may assigned as shown below:

.. code-block:: python

   # Assign the 2Pop model, from model file
   model_2pop =  model_file[‘2Pop’]

A **ModelFile** class object also has two additional attributes:

*ind_file: If created, the filename of a file containing all the unique individuals found within all models stored within the **ModelFile**.
*exclude_file: If created, the filename of a file containing all the unique individuals found within all models stored within the **ModelFile** that do not match a given list of individuals.

These attributes may be populated and their files created using the following functions: **create_ind_file()** and **create_exclude_ind_file()**. The inverse operating can also be completed using **delete_ind_file()** and **delete_exclude_ind_file()**.

***********
Model Class
***********
The **Model** class is used to store model meta-data. The primary attributes of the class are:

*name: The name (as a string) of the model
*tree: The newick tree of the model
*pop_list: A list of the population names within the model
*ind_dict: A dictionary used to store the individuals within each population; the keys of the  dictionary are the population names whereas the values are lists of individuals.
*nind: A dictionary used to store the number of individuals within each population; the keys of the  dictionary are the population names whereas the values are individual counts.
*npop: The number of populations within the model
*inds: A list of all individuals in the model


A Model may be created with all primary attributes populated as shown below:

.. code-block:: python

   # Create the model
   model = pgpipe.model.Model("2Pop")
   
   # Assign the model tree
   model.assign_tree("(A,B);")

   # Assign the populations and their individuals
   model.assign_pop("A", ["Ind1", "Ind2", "Ind3"])
   model.assign_pop("B", ["Ind4", "Ind5", "Ind6"])

A **Model** class object also has two additional attributes:

*pop_files: If created, a list of population filenames. Each population file consist of the individuals found within the population.
*ind_file: If created, the filename of a file containing all the unique individuals found within the model.

These attributes may be populated and their files created using the following functions: **create_ind_file()** and **create_pop_files()**. The inverse operating can also be completed using **delete_ind_file()** and **delete_pop_files()**.

Lastly, a **Model** class object masy be assigned to a **ModelFile** class object as shown below:

.. code-block:: python

   # Create ModelFile object
   models = pgpipe.model.ModelFile()

   # Save the model
   models[str(model.name)] = model