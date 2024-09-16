Databases
=========

The pairinteraction software constructs Hamiltonians from matrix elements that are stored in databases. This design enables the inclusion of new atomic species and even molecules into the software without requiring modifying the code of the software itself. To leverage the software's capabilities for a new species, users only need to supply an appropriate database.

The currently available databases are hosted on GitHub and are downloaded by pairinteraction when needed. Within our GitHub organization, we manage two database repositories:

* `database-sqdt`_: States and matrix elements calculated by single-channel quantum defect theory.
* `database-mqdt`_: States and matrix elements calculated by multi-channel quantum defect theory, utilizing our tool :doc:`mqdt`.

.. _database-sqdt: https://github.com/pairinteraction/database-sqdt
.. _database-mqdt: https://github.com/pairinteraction/database-mqdt

To date, our databases cover alkali metal atoms and a few alkaline earth metal atoms. Nonetheless, our generic database format is equally suited for Rydberg states of other atoms and even molecules. If you are interested in contributing to the databases, please contact us.

The Database Format
-------------------

The databases store tables of states and matrix elements of a set of fundamental operators. We have chosen the set of operators such that it allows to efficiently construct more complex operators such as the operator for dipole-dipole interaction.

.. raw:: html

   <p><details>
   <summary><a>SQL File Defining the Tables of the Database</a></summary>

.. literalinclude:: database.sql
   :language: sql

.. raw:: html

   </details></p>

To enhance performance and minimize memory usage, the tables are converted from SQL to the Apache Parquet format prior to being uploaded to GitHub.

Creating Your Own Databases
---------------------------

You can use our database repositories as templates to create and host your own databases. To tell the pairinteraction software the URL to your database repository, adapt the ``database.json`` file. This file is located in the operating system-specific directory for software configurations, e.g., ``~/.config/pairinteraction/database.json`` on GNU/Linux.
