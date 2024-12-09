.. _advanced:

Advanced Installation
=====================

Advanced users, especially those who want to :ref:`contribute to the development of pairinteraction <repository>`, can build the software from source. The source code is available on :github:`GitHub <>` and can be cloned from the pairinteraction repository using the following `git`_ command:

.. code-block:: bash

    git clone --single-branch https://github.com/pairinteraction/pairinteraction.git

The ``--single-branch`` flag is not essential but will speed up the download significantly by omitting all other branches except master.

For building from source we offer two different setups:
    :ref:`Automatic build <automatic>` - This build method builds the C++ backend automatically, while allowing to build the python dependencies manually. Recommended for adjusting only the python code.

    :ref:`Manual build <manual>` - This build method allows to manually build both the C++ and the python backend. It allows to compile specific parts of the code only, define test targets, and others. Recommend for adjusting both the C++ and the python code.

As building from source depends on your operating system, you will need to prepare different tools for your development environment, including a python package manager.

You will find detailed instructions about how to set up your development environment here:
    :ref:`Python setup <python_setup>` - Get instructions on different python package management systems, and how to setup a virtual environment for development

    :ref:`System setup <system_setup>` - Get instructions on the tooling chain to build the C++ backend via the :ref:`automatic <automatic>` or :ref:`manual <manual>` build.

Building from source allows also to build only the documentation. This is especially useful if you only want to document some changes you contributed to the python interface.
You can also find detailed information about :ref:`building the documentation <Documentation>` on the respective subpage.

.. toctree::
   :maxdepth: 2
   :hidden:

   setup/documentation/documentation
   setup/python/python
   setup/system/system
   build/manual/manual
   build/automatic/automatic


.. _git: https://git-scm.com/downloads/
