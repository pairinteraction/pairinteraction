.. _installation:

Installation
============

The pairinteraction software runs on Linux, macOS 13 or later (Intel & Apple Silicon), and Windows. It is compatible
with Python ≥ 3.9.

Installing from PyPI
--------------------

For most users, we recommend installing pairinteraction from the `Python Package Index (PyPI)`_. You can install
pairinteraction from the command line via the pip_ package manager:

.. code-block:: bash

    pip install pairinteraction

.. tip::

    We strongly advise performing the installation inside a virtual environment to avoid conflicts with other Python
    packages. For an easy setup, even on systems where Python is not yet installed, we recommend using uv_. This
    blazingly fast package manager is becoming increasingly popular as an alternative to pip_. You can run the following
    commands to set up uv_ and install pairinteraction in a new virtual environment with a recent version of Python:

    .. tabs::

        .. tab:: macOS and Linux

            .. code-block:: bash

                # install the uv package manager
                curl -LsSf https://astral.sh/uv/install.sh | sh

                # create a new virtual environment in the current directory
                uv venv --python 3.12

                # activate the environment
                source .venv/bin/activate

                # install pairinteraction
                uv pip install pairinteraction

                # deactivate the environment when you are done using pairinteraction
                deactivate

        .. tab:: Windows

            .. code-block:: bash

                # install the uv package manager
                powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"

                # create a new virtual environment in the current directory
                uv venv --python 3.12

                # activate the environment
                .venv\Scripts\activate

                # install pairinteraction
                uv pip install pairinteraction

                # deactivate the environment when you are done using pairinteraction
                deactivate

.. _pip: https://pypi.org/project/pip/

.. _python package index (pypi): https://pypi.org/project/pairinteraction

.. _uv: https://docs.astral.sh/uv/

Next Steps
~~~~~~~~~~

After the installation, we recommend verifying that everything is working:

.. code-block:: bash

    # print the current version of the pairinteraction library
    pairinteraction --version

    # run a module test to check if everything is working correctly
    pairinteraction test

    # download the databases for the species you want to use (e.g. Rb and Yb174_mqdt)
    pairinteraction download Rb Yb174_mqdt

If the download command fails, you can also download the databases manually from github and place them in the cache
directory. For more details on this, see the :ref:`database <database>` section of the documentation.

**Using pairinteraction as a Python Library:** For usage examples of the Python library, visit the :ref:`tutorial-api`
section of the documentation.

**Using the Graphical User Interface (GUI):** The graphical user interface can be started by executing

.. code-block:: bash

    pairinteraction gui

from the command line. This allows you to do some common and basic calculations without writing any code. For more
information on how to use the graphical user interface, visit the :ref:`tutorial-gui` section of the documentation.

Installing from Flathub (Linux only)
------------------------------------

If you are on Linux and only interested in the graphical user interface, you can install pairinteraction as a Flatpak_
package. For the installation, you have to first install Flatpak (follow the `setup guide`_ for your Linux distribution)
and then run:

.. code-block:: bash

    flatpak remote-add --if-not-exists flathub https://dl.flathub.org/repo/flathub.flatpakrepo
    flatpak install flathub org.pairinteraction.Pairinteraction

Afterwards, you can start the GUI with:

.. code-block:: bash

    flatpak run org.pairinteraction.Pairinteraction

.. _flatpak: https://flathub.org/apps/org.pairinteraction.Pairinteraction

.. _setup guide: https://flathub.org/setup

Building from Source
--------------------

For developers and experienced users who want to adjust the source code to their own needs, or want :ref:`to contribute
to the package <getting_started_as_a_contributor>` we recommend building the software from source. Some instructions on
how to do this can be found in the :ref:`advanced`.

.. toctree::
    :hidden:

    advanced/advanced.rst
