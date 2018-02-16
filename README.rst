Ogo
===
Like the `Ogopogo of the Okanagan`_, osteoporosis is elusive.
Contribute to a future where we can predict osteoporotic fracture.

.. _Ogopogo of the Okanagan: https://youtu.be/AbKw44AmHbY

.. image:: https://readthedocs.org/projects/ogo/badge/?version=latest
    :target: http://ogo.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://circleci.com/gh/Bonelab/Ogo.svg?style=svg
    :target: https://circleci.com/gh/Bonelab/Ogo

Developer Install
=================
.. code:: bash
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
