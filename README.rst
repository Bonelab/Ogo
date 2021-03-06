
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

.. |Azure| image:: https://dev.azure.com/babesler/Bone%20Imaging%20Laboratory/_apis/build/status/Bonelab.Ogo?branchName=master
    :target: https://dev.azure.com/babesler/Bone%20Imaging%20Laboratory/_build/latest?definitionId=14&branchName=master


Developer Install
=================
.. code-block:: bash

    # I would recommend you create a virtual environment
    # ... either using conda
    conda create -n env_name python=2   # or 3
    conda activate env_name

    # ... or using virtualenv
    virtualenv env_name
    . env_name/bin/activate

    # Install in an 'editable' format 
    pip install -e .

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
