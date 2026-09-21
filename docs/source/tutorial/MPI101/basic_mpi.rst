The basic MPI program
=====================

Six routines introduce the basic lifecycle: ``MPI_Init``, ``MPI_Comm_size``,
``MPI_Comm_rank``, ``MPI_Send``, ``MPI_Recv``, and ``MPI_Finalize``.

.. code-block:: c

    #include <stdio.h>
    #include <mpi.h>

    int main(int argc, char **argv)
    {
        MPI_Init(&argc, &argv);
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        double start = MPI_Wtime();

        printf("Hello from rank %d of %d\n", rank, size);

        double elapsed = MPI_Wtime() - start;
        printf("Rank %d: elapsed %f seconds\n", rank, elapsed);
        MPI_Finalize();
        return 0;
    }

``MPI_Init(&argc, &argv)`` initialises the MPI environment; the launcher starts
the processes. Its C declaration is:

.. code-block:: c

    int MPI_Init(int *argc, char ***argv);

Initialise MPI before this tutorial's communication routines. Some queries,
such as ``MPI_Initialized``, are permitted before initialisation.

``MPI_Comm_rank`` and ``MPI_Comm_size`` write results through pointers. A
communicator also isolates message matching: equal tags on different
communicators do not make messages match.

Timing and completion
---------------------

``MPI_Wtime`` provides a local wall-clock reading. Subtract readings on the
same rank to measure an interval; do not assume clocks on different ranks
are synchronised. The maximum rank interval is often useful for assessing
whole-job performance.

Rank output can be interleaved. A barrier ensures group participation at that
point; it does not impose a print order.

Complete outstanding communications before ``MPI_Finalize``. In these
World-model examples all ranks participate in finalisation. An unrecoverable
error can use ``MPI_Abort``; ordinary argument validation should lead all ranks
to the same success or failure path.

Continue to :doc:`monte-carlo` for the first exercise.
