MPI-IO: one shared output file
==============================


Every rank owns part of the solution. MPI-IO lets ranks write their portions
directly into one shared file, avoiding a separate assembly step or gathering
the whole grid into root memory. It can function on supported ordinary
filesystems as well as parallel filesystems; filesystem performance matters
at scale. Follow the training site's current storage guidance.

.. list-table::
   :header-rows: 1

   * - Routine
     - Participation here
     - Purpose
   * - ``MPI_File_open``
     - Collective
     - Open the same file with consistent modes
   * - ``MPI_File_set_size``
     - Collective
     - Set length, including truncation of an older larger file
   * - ``MPI_File_write_at``
     - Independent
     - Write at an explicit offset
   * - ``MPI_File_write_at_all``
     - Collective
     - Alternative write coordinating participation
   * - ``MPI_File_close``
     - Collective
     - Close the file


Do not restrict collective file management to rank 0. ``MPI_MODE_CREATE``
creates a missing file; it does not truncate an existing one. The exercise
therefore sets the size explicitly.

.. code-block:: c

    int MPI_File_write_at(MPI_File fh, MPI_Offset offset, const void *buf,
                          int count, MPI_Datatype datatype, MPI_Status *status);

``buf`` supplies the data to write. ``MPI_File_read_at`` is the corresponding
read routine. The output here is native binary doubles without a header;
MPI-IO does not automatically format numbers as text.

Count and offset use different units
------------------------------------

``count`` measures buffer elements of ``datatype``. File offsets measure
etype units relative to the file view. This example retains the default byte
view, so its offsets are byte positions. Revisit the assumption if changing
the view. See the `MPI file-view definition
<https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node350.htm>`_.

Derive the output layout from ownership
---------------------------------------

Write each global row once: exclude internal ghosts and include both physical
boundaries. With :math:`N` equal to ``mesh_size`` and :math:`b` equal to
``sizeof(double)``, the byte position of global row :math:`g` is

.. math::

    \operatorname{offset}(g)=gNb.

The first owned row on rank :math:`r` is :math:`rn+1`, where :math:`n` is
``int_rows``. Rank 0 also writes row 0; the last rank also writes the top
physical boundary and any remaining interior rows. Use ``MPI_Offset``
arithmetic before multiplying rank, row count, mesh size, and element size.

For a 10 by 10 grid on three ranks, eight interior rows split as 2, 2, and 4.
With 8-byte doubles:

.. list-table::
   :header-rows: 1

   * - Rank
     - Global rows
     - Doubles
     - First byte
     - Byte interval
   * - 0
     - 0--2 (including bottom boundary)
     - 30
     - 0
     - ``[0, 240)``
   * - 1
     - 3--4
     - 20
     - 240
     - ``[240, 400)``
   * - 2
     - 5--9 (including top boundary)
     - 50
     - 400
     - ``[400, 800)``


The intervals have no overlap or gaps. A one-rank run writes all rows including
both physical boundaries. Matching the total file size alone cannot establish
correct row placement; inspect the offsets and values too.

Output size and smaller examples
--------------------------------

Only the final grid is written. Its size is ``mesh_size * mesh_size *
sizeof(double)`` bytes, independent of the iteration count and number of ranks.
On the training platform, each double occupies eight bytes:

.. list-table::
   :header-rows: 1

   * - Mesh
     - Output bytes
     - Approximate size
   * - 300 x 300
     - 720,000
     - 703.1 KiB
   * - 100 x 100
     - 80,000
     - 78.1 KiB
   * - 30 x 30
     - 7,200
     - 7.0 KiB


The 300 by 300 example is under 1 MB. To reduce it, change the first solver
argument and the NumPy reader's ``mesh_size`` to the same smaller value.
For example, ``100 1000 Jacobi`` needs ``mesh_size = 100`` in the reader.
The residual changes with mesh size; the 300-grid checkpoint no longer applies.
The program resizes the file on each run.

.. _exercise-2-3:

Exercise 2.3: Write the middle ranks (optional)
-----------------------------------------------

.. admonition:: STOP HERE -- Exercise 2.3
   :class: important

   **Edit:** :download:`day2/laplace_mpi_io.c <../../../../day2/laplace_mpi_io.c>`

   #. Use at least three ranks (four in the command below) to exercise the middle-rank branch.
   #. Find the final ``#Todo`` in the middle-rank ``else`` branch. Exclude both ghost rows from its count.
   #. Calculate the first owned global row and its byte offset. Write from the first owned local row using ``MPI_File_write_at``.

Save your changes, then run from ``day2/``:

.. code-block:: bash

    make io
    mpiexec -np 4 ./laplace_mpi_io 300 1000 Jacobi

**Checkpoint:** The residual is about ``0.031236`` and the 300 by 300 output holds 90,000 doubles. Complete the inspection below before moving to profiling.

Compare your attempt with :download:`the reference solution <../../../../day2/solution/laplace_mpi_io.c>`.


Read the output
---------------

Run in the same ``day2/`` directory, using the mesh size passed to the solver:

.. code-block:: python

    import numpy as np

    mesh_size = 300
    values = np.fromfile("laplace-soln-whole", dtype=np.float64)
    assert values.size == mesh_size * mesh_size, "Unexpected output size"
    solution = values.reshape(mesh_size, mesh_size)
    assert np.isfinite(solution).all(), "Non-finite solution values"
    print(f"Read {solution.shape} grid; range "
          f"[{solution.min():.6f}, {solution.max():.6f}]")

**Checkpoint:** the reader reports ``(300, 300)`` and no assertion fails.
Check row ownership as well; the regression suite checks numerical values.
If skipping this optional exercise, skip both the MPI run and the reader,
then continue to :doc:`../Topic_six/profiling`.
