Lua Interface
=============
The :code:`iceicle` miniapp exposes functionality via input decks written in Lua to allow for dynamic setup of a variety of problems. 
Check the :code:`examples` directory for some example input files. 
To enable this functionality make sure cmake can access Lua development libraries and turn on the cmake option :code:`ICEICLE_USE_LUA`.
This will download the `sol2 <https://github.com/ThePhD/sol2>`_ library through the :code:`FetchContent` mechanism.

Lua input files should return a single table with all of the configuration options. 
These configuration options are split into modules for each overarching table (see below)
For example - the following input file solves a diffusion equation in 1d:

.. code-block:: lua
   :linenos:

   local fourier_nr = 0.0001

   local nelem_arg = 8;

   return {
      -- specify the number of dimensions (REQUIRED)
      ndim = 1,

      -- create a uniform mesh
      uniform_mesh = {
         nelem = { nelem_arg },
         bounding_box = {
               min = { 0.0 },
               max = { 2 * math.pi }
         },
         -- set boundary conditions
         boundary_conditions = {
               -- the boundary condition types
               -- in order of direction and side
               types = {
                  "dirichlet", -- left side
                  "dirichlet", -- right side
               },

               -- the boundary condition flags
               -- used to identify user defined state
               flags = {
                  0, -- left
                  0, -- right
               },
         },
         geometry_order = 1,
      },

      -- define the finite element domain
      fespace = {
         -- the basis function type (optional: default = lagrange)
         basis = "legendre",

         -- the quadrature type (optional: default = gauss)
         quadrature = "gauss",

         -- the basis function order
         order = 2,
      },

      -- describe the conservation law
      conservation_law = {
         -- the name of the conservation law being solved
         name = "burgers",
         mu = 1.0,

      },

      -- initial condition
      initial_condition = function(x)
         return math.sin(x)
      end,

      -- boundary conditions
      boundary_conditions = {
         dirichlet = {
               0.0,
         },
      },

      -- solver
      solver = {
         type = "rk3-tvd",
         dt = fourier_nr * (2 * math.pi / nelem_arg) ^ 2,
         tfinal = 1.0,
         ivis = 1000
      },

      -- output
      output = {
         writer = "dat"
      }
   }

==============
Dimensionality
==============

:code:`ndim` : specifies the number of dimensions.

This will affect the sizes of many input arrays and internally sets the ``ndim`` template parameter. 
For spacetime simulations this is the number of spatial dimensions +1 for the time dimension.

====
Mesh
====

------------
Uniform Mesh
------------

All simulations need a mesh to define and partition the geometric domain. Currently only uniform hyper-cube mesh generation is supported.

:code:`uniform_mesh` : creates a uniform mesh. There are two sets of supported inputs: bounding box and number of elements in each direction, or a table of nodes for each dimension.

**Required Members**
* Option 1: bounding box

   * :code:`nelem` : the number of elements in each direction

     This is a table of size ``ndim`` and is ordered in axis order (x, y, z, ...)

   * :code:`bounding_box` : the bounding box of the mesh 

      * :code:`min` : the minimal corner of the mesh - if the mesh is the bi-unit hypercube, this is the :math:`\begin{pmatrix} -1 & -1 & ... \end{pmatrix}^T` corner of the mesh

      * :code:`max` : the maximal corner of the mesh - if the mesh is the bi-unit hypercube, this is the :math:`\begin{pmatrix} 1 & 1 & ... \end{pmatrix}^T` corner of the mesh

* Option 2: nodes by dimension 
   * :code:`directional_nodes` a table of tables, there must be ``ndim`` tables inside directional_nodes. Each table contains the nodes in that direction. The mesh nodes will be the cartesian product of these arrays.

* :code:`boundary_conditions` : table to define the boundary condition identifiers for the domain faces 

   The tables in this table follow the face number ordering for Hypercube-type elements. 
   Considering a domain :math:`[-1, 1]^n` - this ordering is first the :math:`-1` side for each dimension in order (x, y, z, ...), 
   then the :math:`+1` sides. 

   In 2D this is: -x (left), -y (bottom), +x (right), +y (top).

   * :code:`types` : the :cpp:enum:`iceicle::BOUNDARY_CONDITIONS` type (check for documentation entries with "Lua name" for the case insensitive name)

   * :code:`flags` : an integer flag to pass to the boundary condition. 
     For Dirichlet and Neumann boundary conditions, this is the (0-indexed) index of the value or function 
     in the table where these boundary conditions are configured later.

**Optional Members**

* ``geometry_order`` : the polynomial order of the basis functions defining the geometry 

   (defaults to 1)

----
gmsh
----

The ``gmsh`` table defines a mesh by reading in a gmsh file 

**Required Table Members**

* ``file`` the relative or absolute path to the ``.msh`` file (using MSH format 4.1)

* ``bc_definitions`` The boundary conditions to prescribe for each physical tag. 
   These definitions are tables that are pairs of boundary condition types (:cpp:enum:`iceicle::BOUNDARY_CONDITIONS`) and integer flags 

.. code-block:: lua 

   bc_definitions = {
      -- 1: initial condition (dirichlet, flag 1)
      { "dirichlet", 1 },

      -- 2: right (extrapolation, flag 0)
      { "extrapolation", 0 },

      -- 3: spacetime-future
      { "spacetime-future", 0 },

      -- 4: left (dirichlet, flag 0)
      { "dirichlet", 0 },
   },

.. code-block:: python
   
   import gmsh;

   # mesh size parameter
   lc = 0.005

   # setup the gmsh model
   gmsh.initialize()
   gmsh.model.add("gmsh_model")

   # =========================
   # = make the bounding box =
   # =========================
   bp1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
   bp2 = gmsh.model.geo.addPoint(1, 0, 0, lc)
   bp3 = gmsh.model.geo.addPoint(1, 1, 0, lc)
   bp4 = gmsh.model.geo.addPoint(0, 1, 0, lc)


   bl1 = gmsh.model.geo.addLine(bp1, bp2);
   bl2 = gmsh.model.geo.addLine(bp2, bp3);
   bl3 = gmsh.model.geo.addLine(bp3, bp4);
   bl4 = gmsh.model.geo.addLine(bp4, bp1);



=======
FESpace
=======

:code:`fespace` : this table defines the finite element space 

**Optional Members**

* ``basis`` the basis function type. Current options are:

   * ``"lagrange"`` Nodal Lagrange basis functions (default)

   * ``"legendre"`` Modal Legendre basis functions 

      uses a Q-type tensor product on hypercube elements

* ``quadrature`` the quadrature function type. Current options are:

   * ``"gauss"`` Gauss Legendre Quadrature (default)

* ``order`` the polynomial order of the basis functions (defaults to 0)

.. note::
   This order must be between 0 and MAX_POLYNOMIAL_ORDER (see ``build_config.hpp``) 
   because the polynomials use integer templates for the order to provide optimization opportunities

================
Conservation Law
================

The ``conservation_law`` module is specific to the Conservation Law Miniapp and specifies the conservation law being solved 
and some physical parameters.

----------------
Required Members
----------------

* ``name`` the name of the conservation law, current options are: 
   * ``burgers`` : The viscous burgers equation 

      :math:`\frac{\partial u}{\partial t} + \frac{\partial (a_j u + \frac{1}{2}b_j u^2)}{\partial x_j} = \mu\frac{\partial^2 u}{\partial x^2}`

   * ``spacetime-burgers`` : The viscous burgers equation in spacetime (see :ref:`Spacetime DG`)

      :math:`\frac{\partial u}{\partial t} + \frac{\partial (a_j u + \frac{1}{2}b_j u^2)}{\partial x_j} = \mu\frac{\partial^2 u}{\partial x^2}`

   * ``euler`` the inviscid Euler equations

     :math:`\frac{\partial \mathbf{U}}{\partial t} + \frac{\partial \mathbf{F}_j(\mathbf{U})}{\partial x_j} = 0`

     :math:`\mathbf{U} = \begin{pmatrix} \rho \\ \rho u_i \\ \rho E \end{pmatrix}, \quad \mathbf{F}_j(\mathbf{U}) = \begin{pmatrix} \rho u_j \\ \rho u_i u_j + p \delta_{ij} \\ u_j(\rho E + p) \end{pmatrix}`

   * ``navier-stokes`` The Navier-Stokes equations
     :math:`\frac{\partial \mathbf{U}}{\partial t} + \frac{\partial \mathbf{F}^c_j(\mathbf{U})}{\partial x_j}
     - \frac{\partial \mathbf{F}^v_j(\mathbf{U})}{\partial x_j} = 0`

     :math:`\mathbf{U} = \begin{pmatrix} \rho \\ \rho u_i \\ \rho E \end{pmatrix}, 
     \quad \mathbf{F}^c_j(\mathbf{U}) = \begin{pmatrix} \rho u_j \\ \rho u_i u_j + p \delta_{ij} \\ u_j(\rho E + p) \end{pmatrix},
     \quad \mathbf{F}^v_j(\mathbf{U}) = \begin{pmatrix} 0 \\ \tau_{ij} \\ u_i\tau_{ij} + q_j \end{pmatrix}`

----------------
Optional Members
----------------

* ``source`` Allows the user to specify a source term :math:`f(x) : \mathbb{R}^d \mapsto \mathbb{R}^m` 
  An example source term is as follows 

.. code-block:: lua

        source = function(x, y)
            local gamma = 1.4;
            local Pr = 0.72;
            local mu = 0.1;

            return
                -4 * math.cos(2 * (x + y)),
                -math.cos(2 * (x + y)) * (14 * gamma - 10)
                - math.sin(4 * (x + y)) * (2 * gamma - 2),
                -math.cos(2 * (x + y)) * (14 * gamma - 10)
                - math.sin(4 * (x + y)) * (2 * gamma - 2),
                -math.cos(2 * (x + y)) * (28 * gamma + 4)
                - 4 * gamma * math.sin(4 * (x + y))
                - 8 * mu * gamma / Pr * math.sin(2 * (x + y));
        end,

----------------------------------
Conservation Law Dependent Members
----------------------------------

Depending on the name of the conservation law. Different configuration options will be available to set up 
the conservation law.

For Burgers Equations:

* ``mu`` The viscosity coefficient for burgers equation

* ``a_adv`` a table the size of the number of **spatial** dimensions for the linear advection term :math:`a_j` in burgers equation

* ``b_adv`` a table the size of the number of **spatial** dimensions for the nonlinear advection term :math:`b_j` in burgers equation

.. note::
   The Spacetime burgers equation will have ``ndim-1`` fields because of the one time dimension.

For Navier-Stokes:

* ``viscosity`` This table specifies the viscosity function 

   * ``name`` specifies the type of viscosity function: either ``constant`` for a constant viscosity 
     or ``sutherlands`` for a dimensionless sutherlands law.

   * ``mu`` if ``constant`` viscosity is specified, this sets the viscosity coefficient

=================
Initial Condition
=================

``initial_condition`` (optional) should be a function that takes an ``ndim`` number of arguments
which represent the coordinates in the physical domain and return the value of the solution(s) at that point.

For single equation conservation laws, this will just be a single value.
For multiple equation conservation laws, this should return a table of values of the required size.

If not specified, this initializes the entire domain to 0.

An example for 2D code: 
  
.. code-block:: lua

   initial_condition = function(x, y)
      return (1 + x) * math.sin(y)
   end

===========================
Boundary Conditions Lua API
===========================

``boundary_conditions`` configures different boundary conditions. 

* ``dirichlet`` set Neumann boundary condition for each flag in order 

   These can be real values, or functions that take ``ndim`` arguments (to represent physical space coordinates) 
   and returns the perscribed value at this location


* ``neumann`` set Neumann boundary condition for each flag in order 

   These can be real values, or functions that take ``ndim`` arguments (to represent physical space coordinates) 
   and returns the perscribed value at this location

* ``extrapolation`` Extrapolate the interior value. 

   Assumes that the exterior state and derivatives are equivalent to the interior state

* ``slip wall`` A frictionless wall, a symmetric boundary condition in the velocity components. Note this isn't strictly equivalent to a symmetric boundary condition because the state gradients are taken to be equivalent to the interior state.

======
Solver
======

``solver`` specifies which solver to use when solving the pde. 
Implicit methods assume no method of lines for time, so are either steady state, or spacetime solutions.

* ``type`` the type of solver to use; current named solvers are:

   * ``"newton"``, ``"newtonls"`` : Newtons method with optional linesearch (Implicit)

   * ``"lm"``, ``"gauss-newton"`` : Regularized Gauss-Newton method :ref:`Gauss-Newton`

   * ``"explicit_euler"`` : Explicit Euler's method 

   * ``"rk3-ssp"``, ``"rk3-tvd"`` : Three stage Runge-Kutta explicit time integration. Strong Stability Preserving (SSP) or Total Variation Diminishing (TVD) versions

* ``ivis`` The visualization (output) is run every ``ivis`` iterations of the solver

* ``idiag`` The diagnostic (output) is run every ``idiag`` iterations of the solver 

* ``verbosity`` general control for the verbosity of the solver when iterating

--------------------------
Explicit Solver Parameters
--------------------------

* ``dt`` set the timestep for explicit schemes

* ``cfl`` the CFL number to use to determine the timestep from the conservation law 

.. note::

   ``dt`` and ``cfl`` are mutually exclusive

* ``tfinal`` the final time value to terminate the explicit method

* ``ntime`` the number of timesteps to take 

.. note::
   ``tfinal`` and ``ntime`` are mutually exclusive


--------------------------
Implicit Solver Parameters
--------------------------

* ``tau_abs`` the absolute tolerance for the residual norm to terminate the solve 

* ``tau_rel`` the relative amount to the initial residual norm at which to terminate the solve

The solve stops when the residual is less than :math:`\tau_{abs} + \tau_{rel} ||r_0||`

* ``kmax`` the maximum number of nonlinear iterations to take 

* ``linesearch`` perform a linesearch along the direction calculated by the implicit solver
   * ``type`` the linesearch type. Currently supported types are ``"wolfe"`` or ``"cubic"`` for cubic interpolation linesearch 
     and ``"corrigan"`` for the linesearch described by Ching et al. 2024 Computer Methods in Applied Mechanics and Engineering

   * ``kmax`` (optional) the maximum number of linesearch iterations (defaults to 5)

   * ``alpha_initial`` (optional) the initial multiplier in the direction (1.0 is the full newton step for Newton solver) (defaults to 1)

   * ``alpha_max`` (optional) the maximum multiplierfor the direction (defaults to 10)

   * ``c1`` and ``c2`` (optional) linesearch coefficients (defaults to 1e-4 and 0.9 respectively)

-----------------------
Gauss-Newton Parameters
-----------------------

* ``lambda_u`` regularization value for pde degrees of freedom -- defaults to :math:`10^{-7}`

* ``lambda_lag`` regularization value for lagrangian regularization -- defualts to :math:`10^{-5}`

* ``lambda_1`` regularization value for curvature penalization (not implemented) -- defaults to :math:`10^{-3}`

* ``lambda_b`` regularization value for geometric degrees of freedom -- defaults to :math:`10^{-2}`

* ``form_subproblem_mat`` set to true if you want to explicitly form the subproblem matrix for ``"guass_newton"`` type

======
Output
======

``output`` specifies the output format for solutions (saves in the ``iceicle_data`` directory)

* ``vtu`` Paraview vtu file (2D and 3D only)

.. warning:: 
   vtu field names may not be correctly put in for NS discretizations. 
   This also only works in serial.
   Use 'vtk' instead.

* ``vtk`` paraview pvtu file (works in parallel)

* ``dat`` Space separated values along the solution in 1D (1D only)

===============
Post Processing
===============
Post processing and error analysis is important for verification and validation. 
This module gives the ability to specify analytical solutions, get domain-wide error values, or plot error fields.

To use this module specify parameters in the ``post`` table.

Specify an analytical or baseline solution with the ``exact_solution`` field.
The value in this field can be a string referring to a precoded exact solution, 
or a function that takes in a parameter for each directional coordinate, and returns a value for single equation systems 
or a table for muti-equation systems.

Specify post processing tasks in a table of strings called ``tasks``

* The ``"plot_exact_projection"`` task forms a L2 projection of ``exact_solution`` and outputs it using the writer specified in the ``output`` table

* The ``"l2_error"`` task performs a L2 error analysis using ``exact_solution``

--------------------------
:math:`L^2` Error Analysis
--------------------------

Given the exact solution :math:`u^{exact} : \mathbb{R}^d \mapsto \mathbb{R}^m`, 
the :math:`L^2` norm of the error :math:`e` is:

.. math:: 

   ||e||_{L^2} = \int_\Omega \sqrt{ (u - u^{exact}(\mathbf{x}))^2 } \;d\mathbf{x} 

For vector valued equations (:math:`m>1`), we compute the error component-wise and then take the vector 2-norm of the error vector.
