Exercise roadmap
================

Complete each stop after reading its associated concept section. File paths
are relative to the repository root. Download links are for reference; run
from the full repository so headers, shared C files, and Makefiles are present.

.. list-table::
   :header-rows: 1

   * - Stop
     - File
     - Task
   * - :ref:`Exercise 1.1 <exercise-1-1>`
     - :download:`day1/MC_pi.c <../../day1/MC_pi.c>`
     - Read and run
   * - :ref:`Exercise 1.2 <exercise-1-2>`
     - :download:`day1/mpi_deadlock.c <../../day1/mpi_deadlock.c>`
     - Run, diagnose, fix, and rerun
   * - :ref:`Exercise 1.3 <exercise-1-3>`
     - :download:`day1/laplace_mpi_blocking.c <../../day1/laplace_mpi_blocking.c>`
     - Complete the missing halo exchange
   * - :ref:`Exercise 1.4 <exercise-1-4>`
     - :download:`day1/laplace_mpi_nonblocking.c <../../day1/laplace_mpi_nonblocking.c>`
     - Combine four requests and wait for completion
   * - :ref:`Exercise 1.5 <exercise-1-5>`
     - :download:`day1/laplace_mpi_persistent.c <../../day1/laplace_mpi_persistent.c>`
     - Initialise, start, and free persistent requests
   * - :ref:`Exercise 2.1 <exercise-2-1>`
     - :download:`day2/laplace_mpi_collective.c <../../day2/laplace_mpi_collective.c>`
     - Gather residual contributions and sum at root
   * - :ref:`Exercise 2.2 <exercise-2-2>`
     - :download:`day2/laplace_mpi_win.c <../../day2/laplace_mpi_win.c>`
     - Combine two ghost rows into one window
   * - :ref:`Exercise 2.3 <exercise-2-3>`
     - :download:`day2/laplace_mpi_io.c <../../day2/laplace_mpi_io.c>`
     - Complete middle-rank output; optional
   * - :ref:`Exercise 2.4 <exercise-2-4>`
     - Four reference files listed in the exercise
     - Run and compare profiles; no C edits

The five Day 1 exercises establish the halo-exchange patterns used on Day 2.
For coding exercises, compare with the same filename in that day's
``solution/`` directory. Those files retain the starter's naming and shared
code while implementing the requested changes.
