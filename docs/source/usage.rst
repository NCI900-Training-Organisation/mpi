Running and checking the labs
=============================

Start with :doc:`start` and follow the :doc:`exercises`. Commands below assume
a complete repository and a compatible compiler/MPI environment. Shared
``mesh.c``, ``solver.c``, and headers need no exercise edits.

Build a completed reference example
-----------------------------------

From the repository root:

.. code-block:: bash

    cd day1/solution
    make blocking
    mpiexec -np 2 ./laplace_mpi_blocking 18 10 Jacobi

Arguments are ``mesh_size max_iter Jacobi``. The first counts points per
dimension including boundaries, and the second is a fixed iteration count,
not a convergence tolerance. The printed scaled residual is
:math:`h^2\|Au-b\|_2` and should be about ``0.488953`` for this small case.
All six reference variants should agree to numerical tolerance.

Common issues
-------------

* **Compilation stops at a TODO:** finish the exercise and remove its marker.
* **Wrong source or missing Makefile:** check the working directory; Day 2
  runs from ``day2/``, not ``day1/``.
* **Undefined notebook variables:** use the setup in :doc:`start` to define
  ``repo_root``, ``Path``, and ``nprocs``.
* **Only sending lines in the deadlock exercise:** they precede the actual
  sends. The timeout ends the hung run; the fix must produce received lines.
* **Unexpected residual:** match grid size, iterations, and method before
  comparing. For ``300 1000 Jacobi`` the expected value is about ``0.031236``.
* **Unexpected binary output shape:** match the NumPy reader's ``mesh_size``
  to the solver argument and read from the same working directory.
* **No mpiP reports:** check the mpiP module and link flags, then inspect the
  newly retained report directory.

Output and profiling
--------------------

``laplace-soln-whole`` holds the final grid as native doubles, with no header.
At size 300 it is 720,000 bytes on the training platform. See
:doc:`tutorial/Topic_five/mpi_io` for row ownership, smaller sizes, and the reader.
See :doc:`tutorial/Topic_six/profiling` for the four comparison runs and worksheet.

Regression checks
-----------------

From the repository root:

.. code-block:: bash

    python3 tests/check_mpi.py

The checks build references in temporary directories and validate residuals,
uneven decompositions, one-rank runs, invalid arguments, and MPI-IO contents
and resizing. They require a working MPI launcher with permission to use
local communication sockets. They do not complete or run the starter TODOs.

Build this documentation
------------------------

From the repository root, in a Python environment with the docs dependencies:

.. code-block:: bash

    python3 -m pip install -r docs/requirements.txt
    sphinx-build -W --keep-going -b html docs/source docs/build/html

Open ``docs/build/html/index.html``. This builds documentation without executing
MPI jobs or notebook cells.
