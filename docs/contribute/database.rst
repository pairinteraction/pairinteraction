.. _database:

Database Format
===============

The PairInteraction software constructs Hamiltonians from matrix elements that are stored in databases. This design
enables the inclusion of new atomic species and even molecules into the software without requiring modification of the
code of the software itself.

The currently available databases are hosted on GitHub and are downloaded by PairInteraction when needed. We manage two
database repositories:

- database-sqdt_: States and matrix elements calculated by single-channel quantum defect theory, utilizing our tool
  ryd-numerov_.
- database-mqdt_: States and matrix elements calculated by multi-channel quantum defect theory, utilizing our tool
  MQDT.jl_.

.. _database-mqdt: https://github.com/pairinteraction/database-mqdt/releases

.. _database-sqdt: https://github.com/pairinteraction/database-sqdt/releases

.. _mqdt.jl: https://github.com/pairinteraction/MQDT.jl

.. _ryd-numerov: https://github.com/pairinteraction/ryd-numerov

You can also manually download the databases from database-sqdt_ and database-mqdt_. Simply download the zip file of the
species you are interested in and extract it. Alternatively, you can create your own custom database (see
https://github.com/pairinteraction/database-sqdt for an example on how to create a database for PairInteraction). Then
move the extracted folder (e.g. ``Rb_v1.2``) to the cache directory of PairInteraction inside the
``pairinteraction/database/tables/`` directory. On Linux, this is usually located at
``~/.cache/pairinteraction/database/tables/``. For other platforms, you can use our command line interface to find the
directory:

.. code-block:: bash

    pairinteraction paths

Your final directory structure should look like this:

.. code-block:: bash

    ~/.cache/pairinteraction/database/tables/
    ├── misc_v1.2
    │   └── wigner.parquet
    ├── Rb_v1.2
    │   ├── states.parquet
    │   ├── matrix_elements_d.parquet
    │   ├── matrix_elements_...
    ├── ...

Note that in addition to the species-specific databases, you also always need to download the ``misc`` folder, which
contains a database with Wigner symbols.

The databases store tables of states and matrix elements of a set of fundamental operators. We have chosen the set of
operators such that it allows to efficiently construct more complex operators such as the operator for dipole-dipole
interaction.

.. raw:: html

    <p><details>
    <summary><a>SQL File Defining the Tables of the Database</a></summary>

.. literalinclude:: ../../src/pairinteraction/database/database.sql
    :language: sql

.. raw:: html

    </details></p>

To enhance performance and minimize memory usage, the tables are converted from SQL to the Apache Parquet format prior
to being uploaded to GitHub.
