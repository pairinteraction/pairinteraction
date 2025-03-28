Style Guide
===========

This guide outlines the coding style and conventions used in the project.

We highly recommend using ruff_ for formatting Python code and clang-format_ for C++ code to ensure consistency. The
easiest way to use these tools is to integrate them as pre-commit hooks into your local repository. To do this, install
the clang-format_ tool and the python package pre-commit_ and execute the following command in the root directory of
your repository:

.. code-block:: bash

    pre-commit install

This adds a pre-commit hook that automatically formats your code and conducts style checks before each commit (if you
are interested in the set of rules that we are using for formatting the code, take a look at the ``[tool.ruff]`` section
in the :github:`pyproject.toml <tree/master/pyproject.toml>` and the :github:`.clang-format <tree/master/.clang-format>`
file). For manual checks at any time, execute:

.. code-block:: bash

    pre-commit run --all-files

.. _clang-format: https://clang.llvm.org/docs/ClangFormat.html

.. _pre-commit: https://pre-commit.com

.. _ruff: https://docs.astral.sh/ruff/

While this fixes a lot of formatting issues automatically, there are still some rules that need to be followed manually
and are described in the following.

Naming Conventions for Python and C++
-------------------------------------

- Use descriptive names! Explicit is always better than implicit (for example, do do name a Hamiltonian ``h`` but
  ``hamiltonian``). The only allowed exceptions to this rule are the names of loop variables or variables in list
  comprehensions.
- Use ``snake_case`` (letters lowercase and words separated by underscores) for methods, variables and functions.
  Moreover, we use snake_case for the value of literals in Python.
- Use ``CamelCase`` (first letter of each word capitalized) for classes.
- Use ``ALL_CAPS`` (all uppercase with words separated by underscores) for global constant variables, mainly magic
  numbers that are only used internally. Moreover, we use ALL_CAPS for enum values in C++.
- Sometimes, it is not directly clear whether a method should be named in *singular* or *plural* form, for example,
  ``get_overlap()`` vs. ``get_overlaps()``. We strictly stick to the following rule: If a method returns a list or
  matrix of objects (such as a list of overlaps), we use the plural form. If a method returns a single object, we use
  the singular form.

File and Folder Names
---------------------

**C++**

- Use ``snake_case`` for folders. Files are named after the class/function they contain. Thus, we use ``CamelCase`` for
  file names containing classes and ``snake_case`` for file names containing functions. Note that most of the time, each
  class has its own file.
- Separate the code into source files (``src/cpp/src/**/*.cpp``) and header files
  (``src/cpp/include/pairinteraction/**/*.hpp``).

**Python**

- Always use ``snake_case`` for folders and files. It is common to have multiple classes in one file if they are closely
  related.
- The code that wrapps the C++ backend is located in the private module ``src/pairinteraction/_wrapped/`` and made
  available to the user via the public modules ``src/pairinteraction/real.py`` and ``src/pairinteraction/complex.py``.
- If you want to add a new feature purely in Python, add a new file to ``src/pairinteraction/`` (if the feature does not
  require much code) or create a new folder in ``src/pairinteraction/`` with multiple files (if the feature is more
  complex).

Type Annotations
----------------

To avoid bugs, make the code more readable, and help with code completion, we use strict typing. Avoid ``void`` and
``dynamic_casts`` in C++. Use type annotations in Python.
