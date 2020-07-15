========================
Popgen Pipeline Platform
========================

.. only:: not html
 
   ------------
   Introduction
   ------------

The Popgen Pipeline Platform (PPP) is a software platform with the goal of reducing the computational expertise required for conducting population genomic analyses. The PPP was designed as a collection of scripts that facilitate common population genomic workflows in a consistent and standardized environment. Functions were developed to encompass entire workflows, including: input preparation, file format conversion, various population genomic analyses, output generation, and visualization. By facilitating entire workflows, the PPP offers several benefits to prospective end users - it reduces the need of redundant in-house software and scripts that would require development time and may be error-prone, or incorrect, depending on the expertise of the investigator. The platform has also been developed with reproducibility and extensibility of analyses in mind.

The PPP was written using the Python programming language and designed to operate using either Python 2.7 or 3.7. However, as `Python 2 will no longer be maintained past January 1, 2020 <https://www.python.org/dev/peps/pep-0373/>`_ we strongly recommend using Python 3. We designed the PPP as a collection of modular functions that users may combine to generate a wide variety of analyses and pipelines. The functions within the PPP are also seperated into four groups: core VCF-based functions; optional BED/STAT file functions; file conversion functions; and analysis functions (Figure 1). 

.. image:: PPP_assets/PPP_Pipeline_Figure.png
   :scale: 50 %
   :align: center

.. centered::
   Figure 1: Structure of the PPP

The core functions of the PPP were designed to operate using VCF-based files primarily due to frequent support for the format among publicly available datasets and population genomics software. Most users will begin their pipelines with these core functions before moving onto an analysis function. Please note that most analysis functions require a preceding file conversion function to operate. 

Please Note: This documentation is currently being devloped and will be updated freqeuntly in the coming days

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   :hidden:

   PPP_pages/install
   PPP_pages/examples
   PPP_pages/functions
   PPP_pages/input_file_generators
   PPP_pages/analyses
   PPP_pages/utilities
   PPP_pages/model
   PPP_pages/contact
   PPP_pages/citations