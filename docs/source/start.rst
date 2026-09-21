.. _getting-started:

Getting started
===============

Open the teaching editions in the complete repository:

* ``day1/mpi_lab1_teaching.ipynb`` -- foundations and point-to-point communication.
* ``day2/mpi_lab2_teaching.ipynb`` -- collectives, RMA, MPI-IO, and profiling.

Download copies for reference:

* :download:`Day 1 teaching notebook <../../day1/mpi_lab1_teaching.ipynb>`
* :download:`Day 2 teaching notebook <../../day2/mpi_lab2_teaching.ipynb>`

Keep the notebooks in their respective day directories. A notebook download
alone does not contain the C exercises, Makefiles, or shared solver files;
use the full training repository to run the labs. The original
``mpi_lab1.ipynb`` and ``mpi_lab2.ipynb`` remain available as shorter lab editions.

Environment
-----------

You need a C compiler, ``make``, an MPI implementation providing ``mpicc`` and
``mpiexec``, and JupyterLab with a Python kernel. Day 2 also uses NumPy; profiling
requires a compatible mpiP installation and its link dependencies.

On Gadi, use a compute allocation or an ARE JupyterLab session. Load the
instructor's compiler, MPI, and mpiP modules before starting Jupyter so the
kernel inherits the environment. Use the same MPI installation to compile
and launch. Choose ranks within your allocated CPUs.

Working directory and notebook variables
----------------------------------------

Adapt site-specific paths in the setup cell to your clone. Day 1 commands run
from ``day1/`` and Day 2 commands from ``day2/``. The following setup also defines
``repo_root``, ``Path``, and ``nprocs``, which the run and profiling cells use:

.. code-block:: python

    import os
    from pathlib import Path

    repo_root = Path("/path/to/mpi").expanduser().resolve()
    day = "day1"  # Change to "day2" in the Day 2 notebook.
    nprocs = 4    # Do not exceed your allocated CPUs.
    os.chdir(repo_root / day)
    print(f"Working directory: {Path.cwd()}; MPI ranks: {nprocs}")

Replace ``/path/to/mpi`` before running. The deadlock exercise always uses two
ranks. Solver exercises require at least two interior rows per rank:
``mesh_size >= 2 * nprocs + 2``. The MPI-IO exercise needs at least three ranks
to execute the middle-rank branch; four ranks work for the supplied examples.

Compiler wrappers
-----------------

MPI wrappers call an underlying compiler and add the MPI include paths and
link options. ``mpicc`` is the C wrapper for the MPI installation loaded in
your environment; it is not necessarily GCC. Intel MPI's ``mpiicc`` defaults
to the classic Intel ``icc`` compiler, while ``mpiicx`` uses ``icx``.
An Open MPI installation built with Intel compilers still uses ``mpicc``.

For this Open MPI training:

.. code-block:: bash

    mpicc --showme:command
    mpicc --showme:compile
    mpicc --showme:link

See the `Open MPI wrapper reference <https://docs.open-mpi.org/en/v5.0.x/man-openmpi/man1/ompi-wrapper-compiler.1.html>`_
and `Intel MPI compiler commands <https://www.intel.com/content/www/us/en/docs/mpi-library/developer-reference-linux/2021-16/compiler-commands.html>`_.

How to use an exercise stop
---------------------------

#. Read the concept section and pause at its numbered exercise.
#. Open the named C file in JupyterLab's text editor.
#. Replace the intentional ``#TODO`` markers as you complete the task.
#. Save the file, then run the notebook cell or documented terminal command.
#. Check the result before continuing; consult the solution after your attempt.

The deadlock exercise has no TODO marker: run it unchanged before fixing it.
``mesh.c``, ``solver.c``, and their headers are supplied infrastructure.
