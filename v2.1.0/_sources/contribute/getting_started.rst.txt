.. _getting_started_as_a_contributor:

Getting Started as a Contributor
================================

The pairinteraction software has greatly benefited from the contributions of our community. If you're also interested in
enhancing pairinteraction, this guide is for you. Don't hesitate to reach out to the maintainers for any questions or
further assistance.

Ways to Contribute
------------------

Reporting Issues
    Encountered a bug or have a suggestion? Help us improve by submitting an issue on :github:`our GitHub issue page
    <issues>`. Your input is invaluable for ongoing improvements and bug fixes.

Contributing New Quantum Defects
    Precise quantum defects are crucial for our software's accuracy. We particularly seek improvements for species with
    two-valence electrons, where data remains sparse. If you've conducted research that provides new quantum defects,
    sharing your data would significantly benefit the community. To ensure that your research receive appropriate
    recognition, our :doc:`README <../index>` promotes the citation of publications that have contributed quantum
    defects.

Writing Tutorials
    Share your knowledge by writing tutorials. If you've published research utilizing pairinteraction, consider
    contributing a tutorial. This not only aids others in replicating your results but also enhances the visibility of
    your work. Your tutorial can serve as a practical guide for applying pairinteraction to solve complex problems.

Developing Features or Resolving Bugs
    Community contributions in the form of new features or bug fixes are highly appreciated. They not only alleviate the
    workload on our maintainers but also enrich the software with diverse expertise.

    You can either implement new features in the C++ backend and write Python bindings and wrappers for them (as an
    example, see the :github:`function for diagonalizing the Hamiltonians of multiple systems in parallel
    <tree/master/src/cpp/src/diagonalize/diagonalize.cpp>`), or write them purely in Python (as an example, see the
    :github:`perturbative calculation of dispersion coefficients and effective Hamiltonians
    <tree/master/src/pairinteraction/perturbative>`).

Making Changes to the Repository
--------------------------------

If you're interested in creating tutorials, adding features, or addressing bugs, and want to make changes to the
pairinteraction repository, this guide is for you. It assumes a basic familiarity with Git and GitHub. For newcomers, we
encourage exploring `GitHub's educational resources`_.

.. _github's educational resources: https://docs.github.com/en/get-started

1. Fork and Clone the Repository
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Start by forking the pairinteraction repository to your GitHub account by clicking the :github:`"Fork" <fork>` button on
the repository's page. This action creates your own version of the repository, enabling you to perform changes freely.

Once forked, clone the repository to your local machine to start working on the files directly.

.. code-block:: bash

    git clone https://github.com/YourGitHubUsername/pairinteraction.git

Familiarize yourself with the repository's architecture. The software is divided into a :github:`C++ backend
<tree/master/src/cpp>` with Python bindings and a :github:`Python library <tree/master/src/pairinteraction>`. The Python
library makes the Python bindings accessible and add additional functionality like the perturbative calculation of
effective Hamiltonians.

2. Set Up a Development Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before implementing any changes, :ref:`set up a development environment and build the software from source <advanced>`
to guarantee the build system is functioning on your computer. In the following we provide a brief summary of the
necessary steps:

- Install the :ref:`tools <system_setup>` for running the build system and a package manager for managing the Python
  dependencies. In addition, you need to install the dependencies of the C++ backend:

  If you are using GNU/Linux, instructions for installing dependencies with ``apt`` / ``yum`` / ``zypper`` can be found
  in the Dockerfiles that are located in the :github:`docker branch <tree/docker/docker>`.

  If you are using OS X, take a look at the ``brew install`` instructions in :github:`our GitHub workflow
  <tree/master/.github/workflows/cpp-backend.yml>`.

  If you are using Windows, you can use ``vcpkg`` to install the dependencies. To set up vcpkg execute the following
  commands in the powershell:

  .. code-block:: bash

      git clone https://github.com/microsoft/vcpkg.git C:\path\to\vcpkg
      C:\path\to\vcpkg\bootstrap-vcpkg.bat
      $env:VCPKG_ROOT = "C:\path\to\vcpkg"
      $env:PATH = "$env:VCPKG_ROOT;$env:PATH"
      $env:CMAKE_TOOLCHAIN_FILE = "$env:VCPKG_ROOT\scripts\buildsystems\vcpkg.cmake"

  Then, to install the dependencies of the C++ backend from :github:`our configuration file <tree/master/vcpkg.json>`,
  run the following command in the root directory of the pairinteraction repository:

  .. code-block:: bash

      vcpkg install --triplet x64-windows

- After installing the dependencies and activating a python environment, you have two options to build the software. You
  can either build the complete software using ``pip``:

  .. code-block:: bash

      pip install -e .[tests,docs]

  Or, you can build solely the C++ backend using ``cmake``. This manual approach is recommended if you are planning to
  contribute to the C++ backend because it allows for a faster build and more fine-grained control.

  .. code-block:: bash

      pip install -r .build_requirements.txt
      mkdir build
      cd build
      cmake ..
      cmake --build .

To ensure your code adheres to the project's coding standards, we highly recommend using the `pre-commit tool`_. Once
you've installed this tool, integrate it as pre-commit hook into your local repository with the following command:

.. code-block:: bash

    pre-commit install

This automatically formats your code and conducts style checks before each commit. For manual checks at any time,
execute:

.. code-block:: bash

    pre-commit run --all-files

.. _pre-commit tool: https://pre-commit.com

3. Implement, Test, and Document Your Changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After applying your changes, run our test cases to ensure that the software is still working. If you built the software
via ``pip``, run the following command (the virtual environment must be activated):

.. code-block:: bash

    pytest

If you used ``cmake``, execute the command below in your build directory to run all C++ tests:

.. code-block:: bash

    ctest -C RelWithDebInfo

If you added new features, consider writing tests to validate their functionality and a tutorial to demonstrate their
usage.

4. Commit and Push
~~~~~~~~~~~~~~~~~~

With successful testing and having added some documentation, commit your changes and push them to your fork (if you are
working on multiple different features, consider creating a new branch for each feature; otherwise, you can commit
directly to the master branch of your fork).

.. code-block:: bash

    git add Path/To/ModifiedFiles
    git commit -m "Your commit message"
    git push

5. Submit a Pull Request
~~~~~~~~~~~~~~~~~~~~~~~~

Finally, initiate a pull request to merge your contributions with the main repository. From the main repository page, go
to the :github:`"Pull requests" <pull>` page, and click the :github:`"New pull request" <compare>` button to compare
your fork to the original pairinteraction repository. After reviewing your changes, submit the pull request for
approval.
