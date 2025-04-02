Style Guide
===========

This guide outlines the coding style and conventions used in the project.

We highly recommend using ruff_ for formatting Python code and clang-format_ for C++ code to ensure consistency. The
easiest way to use these tools is to integrate them as pre-commit hooks into your local repository. To do this, install
the clang-format_ tool and the python package pre-commit_ and execute the following command in the root directory of
your repository:

.. code-block:: bash

    pre-commit install

This adds a pre-commit hook that automatically formats your code and conducts style checks before each commit. The
pre-commit hook applies most `ruff rules <https://docs.astral.sh/ruff/rules/>`_ (see the ``[tool.ruff]`` section in the
:github:`pyproject.toml<tree/master/pyproject.toml>` file) and some C++ clang rules (see the
:github:`.clang-format<tree/master/.clang-format>` file).

For manual checks at any time, execute:

.. code-block:: bash

    pre-commit run --all-files

.. _clang-format: https://clang.llvm.org/docs/ClangFormat.html

.. _pre-commit: https://pre-commit.com

.. _ruff: https://docs.astral.sh/ruff/

While this fixes a lot of formatting issues automatically and enforces some coding rules, there are still some rules
that need to be followed manually and are described in the following.

Naming Conventions for Python and C++
-------------------------------------

- Use descriptive names! Explicit is always better than implicit (for example, don't name a Hamiltonian ``h`` but
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

C++
~~~

- Use ``snake_case`` for folders. Files are named after the class/function they contain. Thus, we use ``CamelCase`` for
  file names containing classes and ``snake_case`` for file names containing functions. Note that most of the time, each
  class has its own file.
- Separate the code into source files (``src/cpp/src/**/*.cpp``) and header files
  (``src/cpp/include/pairinteraction/**/*.hpp``).

Python
~~~~~~

- Always use ``snake_case`` for folders and files. It is common to have multiple classes in one file if they are closely
  related.
- The code that wraps the C++ backend is located in the private module ``src/pairinteraction/_wrapped/`` and made
  available to the user via the public modules ``src/pairinteraction/real.py`` and ``src/pairinteraction/complex.py``.
- If you want to add a new feature purely in Python, add a new file to ``src/pairinteraction/`` (if the feature does not
  require much code) or create a new folder in ``src/pairinteraction/`` with multiple files (if the feature is more
  complex).

Type Annotations
----------------

To avoid bugs, make the code more readable, and help with code completion, we use strict typing. Avoid ``void`` pointers
and ``dynamic_casts`` in C++. Use type annotations in Python.

Python
~~~~~~

**Argument types** ``def my_func(arg1: ArgumentType1, arg2: ArgumentType2) -> ...:``

- Argument types often don't need to be a specific type like `Duck`, but only have to be some `DuckLikeType` (see also
  the `mypy cheat-sheet <https://mypy.readthedocs.io/en/stable/cheat_sheet_py3.html#standard-duck-types>`_).
- List/Tuples of objects with a fixed length should typically be annotated as ``collections.abc.Sequence[Type]``; for a
  variable length use ``collections.abc.Iterable[Type]`` instead (see also `collections.abc abstract base classes
  <https://docs.python.org/3/library/collections.abc.html#collections-abstract-base-classes>`_).
- Numpy arrays should be annotated as ``"pairinteraction.units.NDArray"`` (which is an alias for
  ``numpy.typing.NDArray[Any]``).

**Return types** ``def my_func(...) -> ReturnType:``

- Return types should be annotated as specific as possible. This allows for better type checking and code completion of
  the returned values.
- Lists should be annotated as ``list[Type]``.
- Tuples should be annotated as ``tuple[Type1, Type2]`` or ``tuple[Type1, ...]`` for a variable length.
- Numpy arrays should be annotated as ``"pairinteraction.units.NDArray"`` (which is an alias for
  ``numpy.typing.NDArray[Any]``).
- Pint objects should be annotated as ``"pairinteraction.units.PintFloat"``, ``"pairinteraction.units.PintArray"`` or
  ``"pairinteraction.units.PintSparse"`` for scalar, dense matrix and sparse matrix respectively. These are aliases for
  ``PlainQuantity[float]``, ``PlainQuantity[NDArray]`` or ``PlainQuantity[csr_matrix]``.

Adding License Information
--------------------------

The pairinteraction software is licensed under LGPL v3. Code files should contain the following license header:

.. code-block::

    SPDX-FileCopyrightText: <year of the creation of the file> Pairinteraction Developers
    SPDX-License-Identifier: LGPL-3.0-or-later

If a third-party library is used *and* binary builds of the pairinteraction software include it, license information for
the library should be added to :github:`LICENSES/THIRD-PARTY-LICENSE-<library>.txt <tree/master/LICENSES>`.
