# Introduction to MPI

Start with [Day 1](day1/mpi_lab1.ipynb), then [Day 2](day2/mpi_lab2.ipynb).

For the training session, use the combined **concepts and practice** editions:

- [Day 1 teaching notebook](day1/mpi_lab1_teaching.ipynb): MPI foundations, communication semantics, diagrams, and the point-to-point exercises.
- [Day 2 teaching notebook](day2/mpi_lab2_teaching.ipynb): collectives, RMA, MPI-IO, profiling, and their exercises.


Follow each notebook's exercise roadmap and pause at **Exercise**. Each stop names the exact file to open, what to change, which cell to run, and what to check before continuing.

- Day 1: Monte Carlo pi, domain decomposition, blocking, nonblocking, and persistent communication.
- Day 2: collectives, remote memory access, MPI-IO, and mpiP profiling.
- Each day's `solution/` directory contains completed reference implementations.
- `docs/` contains the accompanying Sphinx tutorial. The standalone `fd_laplace-*` files are older supplementary examples; the notebooks use the modular `laplace_mpi_*` programs.

## Setup

You need a C compiler, `make`, an MPI implementation providing `mpicc` and `mpiexec`, and JupyterLab with a Python kernel. Day 2's output inspection also uses NumPy. Profiling additionally requires mpiP and its link dependencies.

On Gadi, use a compute allocation or an ARE JupyterLab session. Use the project, storage paths, compiler, MPI, and mpiP modules supplied by your instructor. Load those modules before starting Jupyter so the kernel inherits the environment. Compile and run with the same MPI installation, and choose a rank count within your allocated cores.

