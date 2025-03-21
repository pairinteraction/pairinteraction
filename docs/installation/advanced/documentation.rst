.. _documentation:

Building the Documentation
==========================

For building the documentation, we are using Sphinx_. In order to build the documentation, you first need to set up your
:ref:`development environment <system_setup>` and create a :ref:`python environment <python_setup>`. If not stated
otherwise, all commands described are run from the `docs` folder of the pairinteraction repository.

In order to install all dependencies to smoothly run Sphinx_, you should first run

.. code-block:: bash

    pip install .[doc]

In order to build the documentation once, you can run

.. code-block:: bash

    make html

You can then open the documentation by opening the file ``docs/_build/html/index.html`` with your browser.

Alternatively, you can let the documentation automatically rebuild by running

.. code-block:: bash

    make livehtml

This will start a local web server that serves the documentation at ``http://127.0.0.1:8000`` and automatically rebuilds
the documentation whenever you change the source code or a documentation file.

.. note::

    On **Windows**, you have to replace the command ``make`` by ``call make.bat``

.. _sphinx: https://www.sphinx-doc.org/en/master/index.html
