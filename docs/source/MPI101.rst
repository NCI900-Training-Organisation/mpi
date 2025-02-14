MPI 101
---------------

In this section, we will learn the basics of MPI. Many people get by writing parallel code by just using six MPI calls, and this is what we will show you in this section.

Before we start, let's talk about why MPI is the de-facto standard for parallel programming. 
Consider a computing task that takes one processor :math:`t` time to complete where :math:`t` and :math:`N` are two scalars. Now, it's possible to scale the task in two different ways. 

1. Scale to :math:`N` number of processors and it computes a larger task of size :math:`t \times N` in the same :math:`t` time. This is called **weak scaling**.
2. Scale to :math:`N` number of processors and it computes the same task of size :math:`t` in :math:`t/N` time. This is called **strong scaling**. 

We see in both cases, it requires the ability to scale the number of processors to solve the same problem. This is where MPI comes in - utilise multiple processors via communications.


What MPI does and doesn't do
~~~~~~~~~~~~~~~~~
To achieve the scaling mentioned above, MPI provides the following for the programmer:

.. keypoint::
    MPI provides a standard interface for parallel programming that includes:
        1. a communication libraries;
        2. datatypes;
        3. language bindings with C and Fortran;
        4. parallel I/O;
        5. interface for tool support;
        6. and many more

.. note::
    However, what MPI doesn't do for programmers is:
    **It doesn't write or design the parallel implementation for you.**

.. note::
    The MPI users are responsible for the correctness of the implementation of their parallel algorithms. 


