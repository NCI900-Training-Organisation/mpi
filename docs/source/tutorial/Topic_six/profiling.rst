Profiling MPI communication
===========================


The PMPI interface
------------------

MPI provides matching ``PMPI_`` entry points. A profiler can intercept
``MPI_Send``, time a call to ``PMPI_Send``, record information, and return to
the application. Calling PMPI internally avoids recursively entering its own
wrapper. Profiling adds measurement overhead but preserves the application's
communication algorithm.

.. code-block:: text

    application -> MPI_Send wrapper -> PMPI_Send -> MPI implementation
                     start timer      return
                     record call      return to application

Reading an mpiP report
----------------------

mpiP is a lightweight statistical MPI profiler. Its application interval runs
from the end of ``MPI_Init`` to the start of ``MPI_Finalize``. Interpret task
summaries and call-site information together; accumulated time across ranks
is not elapsed job duration. See the `mpiP guide <https://software.llnl.gov/mpiP/>`_.

.. list-table::
   :header-rows: 1

   * - Report item
     - Question
   * - Application time
     - How long is the measured application interval?
   * - MPI time / percentage
     - How much time is inside MPI routines?
   * - Call count
     - How often does the communication occur?
   * - Bytes and message sizes
     - Are transfers small and frequent or larger and fewer?
   * - Variation across ranks
     - Are workloads or arrivals imbalanced?
   * - Dominant call sites
     - Which source locations account for the costs?


High ``MPI_Waitall`` time may reflect transfers, delayed peers, or little
independent work. High collective time can include waiting for slower ranks.
A smaller MPI percentage alone does not prove a faster program.

Keep workload, rank count, placement, allocation, compiler options, and MPI
implementation consistent. Repeat runs. Comparing blocking and nonblocking
at a fixed rank count is a variant comparison; parallel speedup
:math:`S_p=T_1/T_p` compares the same algorithm on one and :math:`p` ranks.

.. _exercise-2-4:

Exercise 2.4: Compare four communication profiles
-------------------------------------------------

.. admonition:: STOP HERE -- Exercise 2.4
   :class: important

   No C edits are required. Use the four completed reference files below.
   Load the instructor's compatible mpiP environment first. Predict the
   dominant MPI calls, run all four cases, and record your observations.

.. list-table::
   :header-rows: 1

   * - Variant
     - Reference source
     - Make target
   * - Blocking
     - :download:`day1/solution/laplace_mpi_blocking.c <../../../../day1/solution/laplace_mpi_blocking.c>`
     - ``blockingP``
   * - Nonblocking
     - :download:`day1/solution/laplace_mpi_nonblocking.c <../../../../day1/solution/laplace_mpi_nonblocking.c>`
     - ``nonblockingP``
   * - Persistent
     - :download:`day1/solution/laplace_mpi_persistent.c <../../../../day1/solution/laplace_mpi_persistent.c>`
     - ``persistentP``
   * - One-sided
     - :download:`day2/solution/laplace_mpi_win.c <../../../../day2/solution/laplace_mpi_win.c>`
     - ``winP``


The notebook's four profiling cells use a 1000 by 1000 grid and 30 iterations.
Run its setup first so ``repo_root``, ``Path``, and ``nprocs`` are defined;
see :doc:`../../start`. To reproduce the same runs from Python:

.. code-block:: python

    import subprocess
    import tempfile
    from pathlib import Path

    repo_root = Path("/path/to/mpi").expanduser().resolve()  # Replace this path.
    nprocs = 4  # Within the compute allocation.

    def profile_example(day, target, executable):
        source_dir = repo_root / day / "solution"
        subprocess.run(["make", target], cwd=source_dir, check=True)
        report_root = repo_root / "day2" / "profile_reports"
        report_root.mkdir(exist_ok=True)
        report_dir = Path(tempfile.mkdtemp(prefix=target + "-", dir=report_root))
        subprocess.run(
            ["mpiexec", "-np", str(nprocs), str(source_dir / executable),
             "1000", "30", "Jacobi"], cwd=report_dir, check=True)
        reports = sorted(report_dir.glob("*.mpiP"))
        if not reports:
            raise RuntimeError(f"No mpiP report in {report_dir}; check setup.")
        for report in reports:
            print(report.read_text())
        print(f"Reports retained in {report_dir}")

    profile_example("day1", "blockingP", "laplace_mpiP_blocking")
    profile_example("day1", "nonblockingP", "laplace_mpiP_nonblocking")
    profile_example("day1", "persistentP", "laplace_mpiP_persistent")
    profile_example("day2", "winP", "laplace_mpiP_win")

The Makefiles link mpiP and its dependencies through ``CLIBS``. Override those
flags for the site's installation if necessary. Keep each run's report
directory so later runs cannot overwrite its evidence.

Results worksheet
-----------------

.. list-table::
   :header-rows: 1

   * - Variant
     - Ranks
     - Mesh / iterations
     - Elapsed time and measurement
     - Dominant MPI call
     - Report / notes
   * - Blocking
     - 
     - 1000 / 30
     - 
     - 
     - 
   * - Nonblocking
     - 
     - 1000 / 30
     - 
     - 
     - 
   * - Persistent
     - 
     - 1000 / 30
     - 
     - 
     - 
   * - One-sided
     - 
     - 1000 / 30
     - 
     - 
     - 


**Checkpoint:** all four runs produce reports, and you can explain the main
costs and differences. There is no predetermined fastest variant.

Final concept checks
--------------------

#. Why must every rank call gather if only root uses its receive buffer?
#. Why does ``MPI_Put`` returning not permit a target to use the incoming data?
#. How do you prove that file writes neither overlap nor omit rows?
#. Why does lower MPI time not necessarily imply greater parallel speedup?

Use collective participation, RMA completion, row ownership, and consistent
elapsed-time measurements to explain your answers. For correctness checks
across all six reference variants, see :doc:`../../usage`.
