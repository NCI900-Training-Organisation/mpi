Model Problem 
---------------
We conclude learning MPI101 by introducing you a model problem. We will use this model problem to demonstrate a wide range of MPI calls in the rest of the tutorial.

.. admonition:: Model Problem
    :class: hint

    **2D Poisson Equation**

    Consider a two-dimensional Poisson equation with Dirichlet boundary condition over a unit square domain :math:`\Omega = [0,1] \times [0,1]`:

    .. math::
        -\nabla^2 u(x,y) = f(x,y), \quad (x,y) \in \Omega
        u = g \quad \text{on} \partial \Omega

    where :math:`f(x,y)` is a given source term and :math:`u(x,y)` is the solution to the Poisson equation. The Dirichlet boundary condition is given by:

    
    