Shared C Routines
=================

``mesh.c`` and ``solver.c`` are present in each lab and solution directory.
Their declarations are in ``mesh.h`` and ``solver.h``.

``init_mesh``
    Initialise the owned rows, physical boundaries, and ghost rows in buffers
    allocated by the caller. The top rank owns any remaining interior rows.

``Jacobi``
    Compute one Jacobi iteration and copy the updated interior into the
    original mesh buffer. Physical boundary values remain fixed.

``Jacobi_int``, ``Jacobi_top``, ``Jacobi_bottom``
    Compute the interior, top, or bottom owned rows into the update buffer.
    The caller completes the relevant halo communications before updating
    the edge rows, then copies all updated interior values into the mesh.

``local_L2_residual``
    Return the local sum of squares of :math:`h^2(Au-b)`. Sum across ranks and
    take a square root to obtain the printed scaled Euclidean residual.
