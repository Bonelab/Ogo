
Ogo
===
Like the `Ogopogo of the Okanagan`_, osteoporosis is elusive.
Help build a fracture free future.

============= ============
     Docs        Tests    
============= ============
|ReadTheDocs|  |Azure| 
============= ============

.. _Ogopogo of the Okanagan: https://youtu.be/aOgKuMV76KM

.. |ReadTheDocs| image:: https://readthedocs.org/projects/ogo/badge/?version=latest
    :target: http://ogo.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. |Azure| image:: https://github.com/Bonelab/Ogo/actions/workflows/main.yml/badge.svg
    :target: https://github.com/Bonelab/Ogo/actions/workflows/main.yml/badge.svg
    :alt: Tests Status


Install
=======

.. code-block:: bash

    # I would recommend you create a virtual environment
    #  using conda
	conda create -n ogo -c numerics88 -c simpleitk -c conda-forge n88tools simpleitk dicom2nifti pbr nose six python=3
    
	# Activate your conda environment, make sure you have conda installed on your system.
    # Best way to install conda is through anaconda
    conda activate ogo
    
    # Install in an 'editable' format 
    pip install -e .

    # You can also run the full install
    python setup.py install


Utilities
=========

ogo.util
----------

:code:`ogo.util.echo_arguments`

- utility function to use with a command-line application, or a function, to automatically print the arguments to the terminal

:code:`ogo.util.Helper`

- numerous helper functions for tasks in the command line functions (many of these are legacy and not used anymore)

:code:`ogo.util.write_csv`

- utility function for writing data stored in a dictionary to a csv file

:code:`ogo.util.graphcuts`

- place holder for installing software for ogoGraphCuts. Read OGO_GRAPHCUTS_INSTALL.txt

ogo.dat
----------

:code:`ogo.dat.OgoMasterLabels`

- a dictionary of labels used by many of the command line interface functions. Any new labels must be defined here. Use ogoPrintLabels to output the dictionary.

:code:`ogo.dat.MassAttenuationTables`

- the NIST table of mass attenuation versus energy levels used for internal calibration.

:code:`ogo.dat.RT_FEMUR_SIDEWAYS_FALL_REF`

- the reference polygon data for alignment prior to finite element analysis.

:code:`ogo.dat.LT_FEMUR_SIDEWAYS_FALL_REF`

- the reference polygon data for alignment prior to finite element analysis.

:code:`ogo.dat.L4_BODY_SPINE_COMPRESSION_REF`

- the reference polygon data for alignment prior to finite element analysis.

ogo.calib
----------

The calibration functionality is templated as described below and is best accessed using CLI :code:`ogoImageCalibration`

:code:`ogo.calib.calibration`

- top class for all types of calibration

:code:`ogo.calib.standard_calibration`

- inherits from :code:`calibration` and implements basic commands

:code:`ogo.calib.internal_calibration`

- function for performing internal calibration; inherits from :code:`standard_calibration` and :code:`calibration`

:code:`ogo.calib.mindways_calibration`

- function for performing phantom calibration; inherits from :code:`standard_calibration` and :code:`calibration`

Command Line Apps
=================

Here is a list of all of the command-line apps that get installed along with the ogo package.
For detailed usage instructions, type the command followed by :code:`-h` into the terminal
with the :code:`ogo` environment activated.

.. list-table::
   :widths: 25 100
   :header-rows: 1

   * - Command
     - Description
   * - :code:`ogoAnalyzeBMD`
     - measure bone mineral density for each label provided from a calibrated CT scan
   * - :code:`ogoCreateDatasetJSON`
     - generates JSON for defining raw CT, labelled CT, and unlabelled CT for machine learning (nnUNet)
   * - :code:`ogoGenerateFEM`
     - prepares a label (femur, spine) for finite element analysis using FAIM (https://bonelab.github.io/n88/)
   * - :code:`ogoGraphCuts`
     - wrapper for using GraphCuts and requires user to install compiled software (see ogo.util)
   * - :code:`ogoImageCalibration`
     - perform internal calibration of phantom calibration
   * - :code:`ogoImageCrop`
     - utility to crop an image
   * - :code:`ogoImageExam`
     - utility to examine the header, histogram and dimensions of a NIfTI image
   * - :code:`ogoIsotropicResampling`
     - resample a 3D computed tomography datasets to new dimensions
   * - :code:`ogoMergeLabels`
     - combine labels from multiple images into a single image; useful for working with TotalSegmentator
   * - :code:`ogoPrintLabels`
     - output the labels used for bones and other tissues
   * - :code:`ogoProcrustes`
     - determine whether two datasets have similar bone anatomy as per Procrustes
   * - :code:`ogoReadPickle`
     - read pickle files used in machine learning
   * - :code:`ogoRepairNifti`
     - fix instances of corrupted NIfTI files
   * - :code:`ogoReplaceLabels`
     - replace labels in an image with a new label or erase a label
   * - :code:`ogoValidate`
     - can validate accuracy of labels from machine learning and has basic repair functions
   * - :code:`ogoVisualize`
     - allow either offscreen or interactive visualization of labels
   * - :code:`ogodcm2nii`
     - convert DICOM files to NIfTI images

Developer Notes
===============

Style Guide Enforcement
-----------------------
`flake8` is used for style guide enforcement. You will not be able to merge without passing the style guide.

.. code-block:: bash

    cd ogo
    flake8

Running Tests
-------------
`nose` is used for running tests. You will not be able to merge without your tests passing. And please, do write tests.

.. code-block:: bash

    cd Ogo  # From root directory, not ogo
    nosetests

Building Docs Locally
---------------------
Use `sphinx-build`. This should rather fast.

.. code-block:: bash

    cd Ogo  # From root directory, not ogo
    sphinx-build docs/ docs/_build/html/
