# Slide review and notebook integration

Source: `MPI-NCI-2026.tex`, reviewed on 18 September 2026.

The deck provides a useful progression from MPI process and message semantics to halo exchanges, collectives, RMA, profiling, and file output. The teaching notebooks now put that context immediately before the corresponding exercise, rather than requiring attendees to switch between slides and a lab.

## Deliverables

- [Day 1: Concepts and Practice](../../day1/mpi_lab1_teaching.ipynb)
- [Day 2: Concepts and Practice](../../day2/mpi_lab2_teaching.ipynb)

Both are copies of the existing lab notebooks. Original executable exercise cells and their relative order are retained; Day 2 adds a small Python illustration of floating-point summation order. Original notebooks, C files, and slide source were not edited in this integration task.

The copies contain objectives, a linked contents list, conceptual explanations, embedded diagrams, questions with expandable answers, and exercise-specific guidance. Existing C exercises and solutions remain linked in the same directory. Diagrams are notebook attachments, so concept illustrations do not depend on the slide asset directory. The C source files are still required to do the exercises.

Each notebook now opens with an exercise roadmap. Numbered exercise sections mark the slide exercise stops and name the exact file to open, the TODO to find where applicable, the task, and the expected result. Coding exercises lead directly into their build/run cells, followed by checkpoints linking to the next topic. Reference solutions are in expandable panels. The Monte Carlo reading exercise, deadlock exercise, optional MPI-IO task, and profiling activity explicitly state whether participants should read, edit, or run files.

## Coverage map

| Slide topics | Notebook placement |
|---|---|
| Workshop scope; what MPI is and is not; SPMD; strong and weak scaling | Day 1 introduction and MPI foundations |
| MPI implementations; build environment | Day 1 and Day 2 setup, using instructor-selected module versions |
| Monte Carlo example; basic program; ranks and communicators | Day 1 before the first program |
| Poisson model, stencil, decomposition, owned and ghost nodes | Day 1 model problem and ghost-row table |
| Operations, procedures, C bindings, four stages | Day 1 before blocking communication |
| Message data and envelope; datatypes and matching | Day 1 semantics section |
| Standard, buffered, synchronous, ready sends; deadlock | Day 1 blocking exercise introduction |
| Requests, wait/test, receive completion, overlap | Day 1 nonblocking exercise introduction |
| Persistent initialisation, activation, completion, release | Day 1 persistent exercise introduction |
| Broadcast, gather/scatter, reduce/allreduce, scan, non-associativity | Day 2 collective exercise introduction |
| Origin/target, windows, displacement units, put/get, fence epochs | Day 2 one-sided exercise introduction |
| PMPI interception and mpiP | Day 2 profiling introduction and results worksheet |
| File management, independent/collective writes, offsets, output ownership | Day 2 MPI-IO introduction and worked output layout |

The existing lab order is retained: Day 2 covers MPI-IO before profiling. Both are independently navigable, and MPI-IO is marked optional when time is limited. The detailed version timeline, historical survey statistics, and unstructured-grid showcase are not needed to complete the exercises and were not reproduced. The acknowledgement and author attribution are retained in Day 1.

## Corrections applied in the teaching copies

Line numbers below refer to the reviewed TeX source. The slide source itself is unchanged.

| Slide location | Issue | Treatment in the notebooks |
|---|---|---|
| Parallel Monte Carlo, around line 251 | Resetting the seed inside the loop repeats the same samples; the snippet also lacks terminators | Explain seeding once before sampling and reducing once afterwards; retain the working `MC_pi.c` exercise |
| “Start” MPI process, line 283 | The displayed `MPI_Init` declaration is not valid C | Use `int MPI_Init(int *argc, char ***argv);` and show a complete minimal program |
| MPI initialisation / communicator, lines 291–301 | `MPI_WORLD_COMM` is not the predefined communicator name; the pre-init restriction is too broad | Use `MPI_COMM_WORLD` and distinguish communication calls from permitted pre-init queries |
| Program timing / finalisation | Rank output order, local timing, and collective exit behaviour need qualification | Explain that a barrier does not order prints, measure elapsed intervals on each rank, and exit consistently |
| Discretisation | The forcing term is written as `f(u_ij)` | Retain the corrected lab's spatial forcing and grid-size convention |
| Nonblocking completion, line 1058 | `MPI_Test` omits its `flag` argument | Show the correct C declaration and distinguish the return code from the completion flag |
| Overlap / persistent examples, around lines 1120–1300 | Incomplete request initialisation, syntax errors, and edge computation before a successful completion check | Explain both send and receive lifetimes, require completion before reading ghosts or reusing outgoing rows, and retain the working lab exercise |
| Persistent resource release | An array name passed to `MPI_Request_free` releases only its first handle | Explicitly require freeing every request after its final completion |
| Collective categories and broadcast | “All to all” is used loosely; broadcast's type/count rule is oversimplified | Separate all-rank participation from `MPI_Alltoall`, describe matching type signatures, and state same-order participation rules |
| Gather arguments / residual example | Duplicated argument names and inconsistent residual routine spelling/pointer use | Give the correct gather signature, per-rank receive count, and residual sum-of-squares formula |
| RMA fence, line 1677 | Fences are described as switching a process between origin and target | Explain roles per transfer and access/exposure epochs; a rank can be both within an epoch |
| mpiP, line 1767 | A specific old Gadi version is presented as the available version | Use the instructor's compatible environment and explain the report independently of a pinned cluster version |
| Profiling exercise | Lower MPI time is conflated with speedup | Distinguish whole-program time, accumulated rank time, variant comparisons, and strong scaling |
| Parallel I/O, line 1814 | MPI-IO is presented as impossible without a parallel filesystem | Separate API correctness from scalable filesystem performance |
| File write, line 1899 | The write slide shows `MPI_File_read_at` and labels the write buffer as output | Show `MPI_File_write_at` with a `const void *` input buffer |
| File offsets and format | Offset units and binary formatting need qualification | Explain etypes, the default byte view, application-formatted text, and native-double output |
| File output example | Reusing a smaller output file and single-rank boundary ownership are not covered | Explain collective resizing, exact row ownership, and both physical boundaries in a single-rank run |

## Figures and source completeness

At the initial review, only the TeX file was present in `docs/slides/`, so the teaching copies used figures from `docs/source/figures/` and text diagrams for the RMA and MPI-IO explanations. The full slide build stack has since been supplied. The Day 1 blocking communication section now embeds `docs/slides/Blocking_Send.png` directly; the other illustrations retain their existing sources. No external slide viewer is required for these concepts.

The C prototypes and semantic clarifications were checked against the MPI Forum documentation; profiling terminology was checked against LLNL's mpiP guide. Relevant references are linked next to the corresponding explanations in the notebooks.

## Validation scope

Both copies passed notebook-schema validation and HTML export. The exported HTML was checked for embedded image data, tables, expandable answers, and valid navigation anchors. Relative file links resolve. All 14 nonempty original code cells are retained verbatim and in order across the two copies, and execution outputs are cleared. The added floating-point demonstration was run and produces `1.0` and `0.0` for the two addition orders. Checksums confirm the original notebooks and slide source were unchanged.

Day 1 contains 39 cells, seven embedded figures, and thirteen rendered tables. Day 2 contains 45 cells, three embedded figures, and nine rendered tables. HTML previews were generated under `/tmp/mpi-teaching-preview/` for inspection; the deliverables are the notebooks.

The exercise-navigation update passed schema validation and HTML export, including checks of file links and navigation anchors. All five Day 1 exercise stops and four Day 2 exercise stops were checked against their associated run cells. Promoting the deadlock activity to an exercise preserved executable cell contents and embedded image attachments.

The Day 1 deadlock activity is now **STOP HERE — Exercise 1.2: Diagnose and fix a deadlock**, with subsequent exercises numbered 1.3–1.5. Its runnable code cell compiles `mpi_deadlock.c` and launches exactly two ranks with a five-second timeout and a two-second forced-cleanup grace period. Participants first observe and explain the timeout, then edit the source to use `MPI_Sendrecv` and rerun the same cell. An expandable reference fix and a checkpoint requiring both received-message lines follow the run cell. Notebook validation, HTML export, and shell syntax checks passed; this notebook-cell addition was not executed under MPI.

The MPI exercise cells were not rerun for this content-only integration: they require completing the intentional TODOs, and profiling requires the training machine's mpiP installation. HTML export used the Python lexer fallback because IPython was not installed in the temporary validation environment; notebook execution semantics are unchanged.

The Day 1 deadlock example, `mpi_deadlock.c`, now uses a fixed 16 MiB message and static buffers, with no size argument or allocation code. Its complete source is included in the teaching notebook, together with a bounded run command and a `MPI_Sendrecv` fix. The simplified source passes `-Wall -Wextra -Werror` syntax checks. Earlier runtime checks of the configurable version confirmed a 16 MiB standard-send deadlock and successful completion with `MPI_Sendrecv` on local MPICH; the simplified version has not been rerun.

The Day 1 message-semantics section includes C signatures and argument tables for `MPI_Send` and `MPI_Recv`, checked against the linked Open MPI 5.0.x references. The tables distinguish send counts from receive capacity, input arguments from output buffers/status, and error-code returns from received counts. Notebook schema and HTML rendering checks passed, including all six send arguments and seven receive arguments.

The nonblocking section adds matching seven-argument tables and C signatures for `MPI_Isend` and `MPI_Irecv`, with examples from the starter halo exchange. Linked Open MPI references support the request-handle and buffer-lifetime explanations. Both tables passed notebook schema and HTML rendering checks.
