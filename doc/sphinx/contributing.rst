.. _Contributing:

Contributer's Guide
===================

The pairinteraction software has greatly benefited from the contributions of our community. If you're also interested in enhancing pairinteraction, this guide is for you. Don't hesitate to reach out to the maintainers for any questions or further assistance.

Ways to Contribute
------------------

Reporting Issues
    Encountered a bug or have a suggestion? Help us improve by submitting an issue on :github:`our GitHub issue page <issues>`. Your input is invaluable for ongoing improvements and bug fixes.

Contributing New Quantum Defects
    Precise quantum defects are crucial for our software's accuracy. We particularly seek improvements for species with two-valence electrons, where data remains sparse. If you've conducted research that provides new quantum defects, sharing your data would significantly benefit the community. We highlight such contributions in our :doc:`README <index>` to ensure both the software and your research receive appropriate recognition and citation.

Writing Tutorials
    Share your knowledge by writing tutorials. If you've published research utilizing pairinteraction, consider contributing a tutorial. This not only aids others in replicating your results but also enhances the visibility of your work. Your tutorial can serve as a practical guide for applying pairinteraction to solve complex problems.

Developing Features or Resolving Bugs
    Community contributions in the form of new features or bug fixes are immensely appreciated. They not only alleviate the workload on our maintainers but also enrich the software with diverse expertise. Your contributions can make a significant difference, especially in areas outside our maintainers' specialties.

Contributing to the Repository
------------------------------

If you're interested in creating tutorials, adding features, or addressing bugs, you'll need to make changes to files within the pairinteraction repository.
Below, we offer guidance to help you navigate this process effectively, assuming a basic familiarity with Git and GitHub. For newcomers, we encourage exploring `GitHub's educational resources`_.

.. _GitHub's educational resources: https://docs.github.com/en/get-started

Step-by-Step Instructions
^^^^^^^^^^^^^^^^^^^^^^^^^

1. **Fork the Repository:** Start by forking the pairinteraction repository to your GitHub account. This action creates your own version of the repository, enabling you to perform changes freely. Click the :github:`"Fork" <fork>` button on the repository's page to begin.

2. **Clone and Set Up the Development Environment:** Once forked, clone the repository to your local machine to start working on the files directly. Before implementing any changes, build the software from source to guarantee the build system is functioning on your computer. Follow the :doc:`installation <installation>` section of the documentation for detailed steps, specifically the instructions for a manual build, to prepare your development environment properly.

3. **Understand the Structure:** Familiarize yourself with the repository's architecture. The software is divided into a :github:`C++ backend <tree/master/pairinteraction_backend>` with Python bindings, a :github:`Python library <tree/master/pairinteraction>` for added functionality, and a :github:`graphical user interface <tree/master/pairinteraction_gui>` leveraging the Python library.

4. **Implement Your Changes:** With your setup ready, proceed to make your contributions. This could involve coding new features, fixing bugs, or refining documentation. To ensure your code adheres to the project's coding standards, we highly recommend using the `pre-commit tool`_. Once you've installed this tool, integrate our pre-commit hook into your local repository with the following command:

   .. code-block:: bash

       pre-commit install

   This automatically formats your code and conducts style checks before each commit. For manual checks at any time, execute:

   .. code-block:: bash

       pre-commit run --all-files

.. _pre-commit tool: https://pre-commit.com

5. **Test Thoroughly:** After applying your changes, test them to ensure the software's integrity remains intact. If you've followed the manual build guide of the :doc:`installation <installation>` section, execute the command below in your build directory to run all tests:

   .. code-block:: bash

       cmake --build . --target test

   If you added new features, consider writing tests to validate their functionality.

6. **Commit and Push:** With successful tests, commit your changes and push them to your fork.

7. **Submit a Pull Request:** Finally, initiate a pull request to merge your contributions with the main repository. From the main repository page, go to the :github:`"Pull requests" <pull>` page, and click the :github:`"New pull request" <compare>` button to compare your fork to the original. After reviewing your changes, submit the pull request for approval.
