Introduction to MPI
===================

This NCI workshop combines MPI concepts with C programming exercises over two
half-days. Work through Day 1, then Day 2. Each exercise stop names the source
file to edit, the task, a run command, and a checkpoint before continuing.

You will learn to:

* Divide work between ranks and match messages safely.
* Recognise deadlock and manage blocking, nonblocking, and persistent exchanges.
* Combine distributed results with collectives.
* Exchange data through RMA windows and write a shared file with MPI-IO.
* Compare communication costs using mpiP.

The :doc:`teaching notebooks <start>` combine the concepts and practice in one
place. These pages provide the same learning sequence as a searchable reference.
The workshop routines largely predate MPI 4.0; the history section includes
updates through MPI 5.0. An implementation's version number is distinct from
the version of the MPI standard it supports.

.. toctree::
   :maxdepth: 2
   :caption: Training contents

   start
   exercises
   tutorial
   usage
   api

Acknowledgement
---------------

The National Computational Infrastructure acknowledges the Ngunnawal and
Ngambri people of the Canberra region and all First Nations Australians on
whose lands this training takes place. We pay our respects to Elders past
and present.

Prepared from Frederick Fung's NCI workshop slides and teaching notebooks.
