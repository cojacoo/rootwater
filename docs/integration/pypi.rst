====================
Python Package Index
====================

PyPI
====

Setup
-----
.. note::

    Before you can publish any package to the Python Package Index,
    you'll have to register here: https://pypi.org/account/register/

In a second step, you will have to configure the ``setuptools`` to be able to
pack your project file. The easiest approach is to add a new hidden file called
``.pypirc`` to your home folder. This file should look like:

.. code-block:: text
    :caption: .pypirc

    [disutils]
    index-servers =
        pypi

    [pypi]
    repository: https://www.python.org/pypi
    username: yourname
    password: yourpassword

Replace ``yourname`` and ``yourpassword`` with the username and password you
chose when registering on PyPI.

.. warning::

    Your password is saved as plain text in this file, so make sure that
    nobody else has access to your home folder. Otherwise refer to the PyPI
    documentation for upload procedures without saved password.


Prepare your project
--------------------

.. todo::

    Write this section

Upload your project
-------------------

.. todo::

    Write this section