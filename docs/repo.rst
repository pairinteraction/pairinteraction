Contributing to the Repository
==============================

If you're interested in creating tutorials, adding features, or addressing bugs, you'll need to make changes to files within the pairinteraction repository.
Below, we offer guidance to help you navigate this process effectively, assuming a basic familiarity with Git and GitHub. For newcomers, we encourage exploring `GitHub's educational resources`_.

.. _GitHub's educational resources: https://docs.github.com/en/get-started

Step-by-Step Instructions
-------------------------

1. **Fork and Clone the Repository:** Start by forking the pairinteraction repository to your GitHub account by clicking the :github:`"Fork" <fork>` button on the repository's page. This action creates your own version of the repository, enabling you to perform changes freely.

Once forked, clone the repository to your local machine to start working on the files directly.

    .. code-block:: bash

        git clone https://github.com/YourGitHubUsername/pairinteraction.git

Familiarize yourself with the repository's architecture. The software is divided into a :github:`C++ backend <tree/master/src/cpp>` with Python bindings and a :github:`Python library <tree/master/src/pairinteraction>`. The Python library makes the Python bindings accessible and add additional functionality like a graphical user interface.


2. **Set Up the Development Environment:** Before implementing any changes, build the software from source to guarantee the build system is functioning on your computer. Follow the :doc:`installation <installation>` section of the documentation for detailed steps. In the following we provide a brief summary of the necessary steps:

If you are planning to only contribute to the python package a basic installation should be sufficient.

    * Install `uv`_ then go inside the repository and run

    .. code-block:: bash

        uv venv --python=3.9 .venv
        source .venv/bin/activate
        uv pip install -e .[test]


If you are planning to contribute to the C++ backend it might be beneficial to do a manual build:

    * Install the build dependencies `CMake`_ for running the build system and `uv`_ for managing the Python dependencies. In addition, you need to install dependencies of the C++ backend. If you are using GNU/Linux or OS X, dependencies can be found in the Dockerfiles that are located in the :github:`docker branch <tree/docker/docker>`. If you are using Windows, you can use `VCPKG`_ with :github:`our configuration file <tree/master/vcpkg.json>` to install the dependencies. Afterwards, create and activate a python environment and install the necessary build dependencies

    .. code-block:: bash

        uv venv --python=3.9 .venv
        source .venv/bin/activate
        uv pip install nanobind numpy typing_extensions

    If you want to use mkl you should also run ``uv pip install mkl mkl-devel``.

    * Build the software using CMake:

    .. code-block:: bash

        mkdir build
        cd build
        cmake ..
        cmake --build .

3. **Implement Your Changes:** With your setup ready, proceed to make your contributions. This could involve coding new features, fixing bugs, or refining documentation. To ensure your code adheres to the project's coding standards, we highly recommend using the `pre-commit tool`_. Once you've installed this tool, integrate our pre-commit hook into your local repository with the following command:

    .. code-block:: bash

        pre-commit install

    This automatically formats your code and conducts style checks before each commit. For manual checks at any time, execute:

    .. code-block:: bash

        pre-commit run --all-files

.. _pre-commit tool: https://pre-commit.com

4. **Test Thoroughly:** After applying your changes, test them to ensure the software's integrity remains intact.

If you used the basic python installation you can run the tests with the following command (the virtual environment must be activated):

    .. code-block:: bash

        pytest

If you've followed the manual build guide, execute the command below in your build directory to run all C++ tests:

    .. code-block:: bash

        cmake --build . --target test

If you added new features, consider writing tests to validate their functionality.

5. **Commit and Push:** With successful tests, commit your changes and push them to your fork:

    .. code-block:: bash

        git add Path/To/ModifiedFiles
        git commit -m "Your commit message"
        git push


6. **Submit a Pull Request:** Finally, initiate a pull request to merge your contributions with the main repository. From the main repository page, go to the :github:`"Pull requests" <pull>` page, and click the :github:`"New pull request" <compare>` button to compare your fork to the original. After reviewing your changes, submit the pull request for approval.

.. _cmake: https://cmake.org
.. _uv: https://pypi.org/project/uv/
.. _VCPKG: https://vcpkg.io
