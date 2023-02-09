
Ogo
===
Like the `Ogopogo of the Okanagan`_, osteoporosis is elusive.
Help build a fracture free future.

============= ============
     Docs        Tests    
============= ============
|ReadTheDocs|  |Azure| 
============= ============

.. _Ogopogo of the Okanagan: https://youtu.be/AbKw44AmHbY

.. |ReadTheDocs| image:: https://readthedocs.org/projects/ogo/badge/?version=latest
    :target: http://ogo.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. |Azure| image:: https://github.com/Bonelab/Ogo/actions/workflows/main.yml/badge.svg
    :target: https://github.com/Bonelab/Ogo/actions/workflows/main.yml/badge.svg
    :alt: Tests Status


Developer Install
=================
.. code-block:: bash

    # I would recommend you create a virtual environment
    #  using conda
	conda create -n ogo -c numerics88 -c simpleitk -c conda-forge n88tools simpleitk dicom2nifti pbr nose six python=3
    
	# Activate your conda environment, make sure you have conda installed on your system.
    # Best way to install conda is through anaconda
    conda activate ogo
    
    # Install in an 'editable' format 
    pip install -e .

    # Install batchgenerators for ogoPrintPickle
    pip install batchgenerators

    # You can also run the full install
    python setup.py install

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
