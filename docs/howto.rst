================
Use the Template
================

.. toctree::
    :maxdepth: 2

.. _howto:

How to use this template
========================

.. important::

    Although this template is meant to be used for publishing your code on
    the Python Package index (PyPI_), you cannot install this template using
    pip. Please follow this description


.. _PyPI: https://pypi.python.org

Download the template
=====================

Download and unzip the template from http://github.com/KIT-HYD/package-template
either by downloading it directly or using the command line:

.. code-block:: bash

    wget https://gitub.com/KIT-HYD/package-template/archive/master.zip
    unzip master.zip


.. important::

    Do not clone the repository as it will download the hidden ``.git``
    folder as well. This will prevent from initializing your own git project
    at a later step.

Initialize your Project
=======================

Before you continue, make sure to signup at GitHub, GitLab, BitBucket or any
other repository that supports git. As far as I know, there is no way to
create a new project locally and just push it to the repository. Therefore
initializing your project need some more steps:

Create a new repository on e.g. GitHub. I will use GitHub for this guide. You
will have to use a unique name for your project on account level. However,
if you want to add your package also to the Python Package Index, you will
need a global unique name. So check the PyPI before you choose a name.

Once your repository is created, clone it to your computer, wherever you want
your files to live. On the command line this can be done by:

.. code-block:: bash

    cd to/your/path
    git clone https://github.com/username/reponame

Of course you need to change username and reponame to your project
information. In case you use Windows, you can also use the Desktop
application for this step: https://desktop.github.com.
In case you want to use PyCharm (https://www.jetbrains.com/pycarm ), you can
create a new Project directly from the IDE. This is my preferred option.

Rename the template
===================

After unzipping the package-template you will find a folder called
``package-template-master`` in that directory. Copy the content of this
folder into your newly created github project directory.
Depending on the repository you chose, this project might already contain a
``README.md``, ``LICENSE`` and definitely a ``.gitignore`` file. Make sure to
replace these files with the versions from ``package-template-master``.

Commit changes
==============

Finally you need to add the new files to the GitHub repository. Only added
files will be versioned. At some point you can then `commit` changes. A
commit can be thought of a local snapshot of your project, which can also be
reverted. Once you created a commit, this commit can then be `pushed` to the
remote repository on GitHub. Then all changes are also available on the
online version.

In case you are using the Desktop application, you can add, commit and push
from within the software.

PyCharm offers a VCS menu, where you can choose to commit changes. The dialog
will also show you all unversioned (not added) files. Instead of just commit
your changes, you can also choose the 'commit and push' option.

From the command line these steps can be achieved like:

Inside the GitHub repository, all files can be added by:

.. code-block:: bash

    git add .

Then commit the changes

.. code-block:: bash

    git commit -m "initial commit"

Where the -m flag will set a commit message. It is highly recommended to
always add a very short commit message, to keep track of the changes.

Finally, push the commit to the master branch.

.. code-block:: bash

    git push

Integrations
============

This far, you jst set up a GitHub (or gitlab, bitbucket) repository and
copied the template files into your project.
The main purpose of the template is to make the integration of your project
into helpful developer tools easier. But you will usually have to sign up for
thrid party services. Luckily you can use your GitHub account for most
services, with the exception of PyPI.

.. seealso::

    Refer to the following pages to get started with the different developer
    tools:

        * :doc:`Python Package Index </integration/pypi>`