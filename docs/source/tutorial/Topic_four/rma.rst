One-sided communication (RMA)
=============================


Remote memory access uses an **origin** process to transfer data to or from a
**target** process's exposed memory. A **window** describes that memory and the
participating group. Each transfer does not need a matching receive, but setup
and synchronisation still coordinate the ranks.

.. list-table::
   :header-rows: 1

   * - Operation
     - Data direction
     - Use
   * - ``MPI_Put``
     - Origin to target window
     - Deliver an owned row
   * - ``MPI_Get``
     - Target window to origin
     - Fetch exposed data
   * - ``MPI_Accumulate``
     - Combine into target window
     - Apply an operation such as sum


Origin and target are roles for a transfer. The same process can play both
roles within an epoch; a fence does not switch one role into the other.

Memory and displacement units
-----------------------------

``MPI_Alloc_mem`` allocates storage. Collective ``MPI_Win_create`` exposes
existing storage; it does not allocate it:

.. code-block:: c

    int MPI_Win_create(void *base, MPI_Aint size, int disp_unit,
                       MPI_Info info, MPI_Comm comm, MPI_Win *win);
    int MPI_Put(const void *origin_addr, int origin_count,
                MPI_Datatype origin_datatype, int target_rank,
                MPI_Aint target_disp, int target_count,
                MPI_Datatype target_datatype, MPI_Win win);

.. list-table::
   :header-rows: 1

   * - Quantity
     - Units and meaning
   * - ``base``
     - Start of this rank's exposed storage
   * - ``size``
     - Total exposed size in bytes
   * - ``disp_unit``
     - Byte scale for target displacements
   * - ``target_disp``
     - Offset in units of the target's disp_unit
   * - ``origin_count``, ``target_count``
     - Elements described by the associated datatypes


The target byte address is ``target_base + target_disp * target_disp_unit``.
For ``disp_unit = sizeof(double)``, displacement 8 selects the ninth double,
not byte 8. Origin and target data signatures must be compatible and the
written range must fit in the target window.

Fences and completion
---------------------

An access epoch is an interval for origins to issue RMA operations; an exposure
epoch is the corresponding interval for targets. This lab uses collective fences:

.. code-block:: text

    fence -> issue puts to neighbours -> fence -> use received window data

``MPI_Put`` returning alone does not make the transferred data ready to use at
the target. The closing fence provides the completion and visibility needed
here. Keep outgoing data unchanged until completion too. Use
``MPI_Win_fence(0, win)`` for the basic pattern; the first argument is an
assertion bitmask, not a rank or success flag.

See the `MPI_Put definition <https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node317.htm>`_
and `fence synchronisation <https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node327.htm>`_.

One window for two ghost rows
-----------------------------

Allocate ``2 * mesh_size`` doubles. Sketch which half receives each incoming
row. Derive every put's displacement from the **target's** layout:

.. code-block:: text

    target window: [one ghost row: mesh_size doubles]
                   [other ghost row: mesh_size doubles]

Keep puts to distinct locations. Do not copy uninitialised window memory into
a physical boundary when a neighbour is absent. A one-rank run has no incoming
halos and must preserve both boundaries. After the final epoch, collectively
free the window, then free the separately allocated memory.

.. _exercise-2-2:

Exercise 2.2: Use one window for both ghost rows
------------------------------------------------

.. admonition:: STOP HERE -- Exercise 2.2
   :class: important

   **Edit:** :download:`day2/laplace_mpi_win.c <../../../../day2/laplace_mpi_win.c>`

   #. Find the one-window TODO and replace the two allocations/windows with one layout.
   #. Update puts, target displacements, fences, and copies into local ghost rows.
   #. Preserve physical boundaries; release the window and its memory after use.

Save your changes, then run from ``day2/``:

.. code-block:: bash

    make win
    mpiexec -np 4 ./laplace_mpi_win 300 1000 Jacobi

**Checkpoint:** The residual is about ``0.031236``. Explain each target displacement in double-sized units. Continue to the optional MPI-IO exercise.

Compare your attempt with :download:`the reference solution <../../../../day2/solution/laplace_mpi_win.c>`.
