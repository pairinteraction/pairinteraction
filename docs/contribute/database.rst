.. _database:

Database Format
===============

The pairinteraction software constructs Hamiltonians from matrix elements that are stored in databases. This design
enables the inclusion of new atomic species and even molecules into the software without requiring modifying the code of
the software itself.

The currently available databases are hosted on GitHub and are downloaded by pairinteraction when needed. We manage two
database repositories:

- database-sqdt_: States and matrix elements calculated by single-channel quantum defect theory, utilizing our tool
  ryd-numerov_.
- database-mqdt_: States and matrix elements calculated by multi-channel quantum defect theory, utilizing our tool
  MQDT.jl_.

.. _database-mqdt: https://github.com/pairinteraction/database-mqdt/releases

.. _database-sqdt: https://github.com/pairinteraction/database-sqdt/releases

.. _mqdt.jl: https://github.com/pairinteraction/MQDT.jl

.. _ryd-numerov: https://github.com/pairinteraction/ryd-numerov

The databases store tables of states and matrix elements of a set of fundamental operators. We have chosen the set of
operators such that it allows to efficiently construct more complex operators such as the operator for dipole-dipole
interaction.

.. raw:: html

    <p><details>
    <summary><a>SQL File Defining the Tables of the Database</a></summary>

.. literalinclude:: ../../src/pairinteraction/_wrapped/database/database.sql
    :language: sql

.. raw:: html

    </details></p>

To enhance performance and minimize memory usage, the tables are converted from SQL to the Apache Parquet format prior
to being uploaded to GitHub.
