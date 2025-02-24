Technical Background
====================
Partial differential equations (PDEs) are expressed on a domain of interest in :math:`d` dimensional space. 
This domain :math:`\Omega \subset \mathbb{R}^d` is then tesselated by non overlapping elements :math:`\mathcal{K}` in a tesselation :math:`\mathcal{T}`.
We also consider the set of element interfaces :math:`\mathcal{E}` such that :math:`\cup_{\Gamma\in\mathcal{E}} = \cup_{\mathcal{K}\in\mathcal{T}} \partial \mathcal{K}`.

The primary goal is to represent the solution to the PDE accurately by a set of functions on the domain. 
We generally consider PDEs in the form of a **conservation law**; many physical processes follow conservation laws.

.. math::
   \nabla\cdot\mathcal{F}(\mathbf{U}) = 0

Where the quantities in :math:`\mathbf{U}` are **conserved quantities**, such as mass, momentum, and energy. 
We denote the number of conserved quantities with :math:`m`. 
Then the **conservative flux** :math:`\mathcal{F} : \mathbb{R}^m \mapsto \mathbb{R}^{m\times d}`.

If the problem is time-dependent, it is common to treat the time dimension separately and use a conservative flux defined on the spatial dimensions :math:`F : \mathbb{R}^m \mapsto \mathbb{R}^{m\times d_x},\;d_x = d - 1`.
In this case our conservation law becomes:

.. math::
   
   \frac{\partial \mathbf{U}}{\partial t} + \frac{\partial F_j(\mathbf{U})}{\partial x_j} = 0

============================
Function Spaces
============================

We want to find functions that best solve the PDE, so we must consider spaces of functions. First we need a few definitions.

A function :math:`f` is :math:`C^0(\Omega)` or **continuous** on a domain :math:`\Omega` if for all :math:`x_0 \in \Omega`

1. :math:`f(x_0)` is defined 
2. :math:`\lim_{x\to x_0}` exists
3. :math:`\lim_{x\to x_0} = f(x_0)`

For more details see the `Mathworld Article <https://mathworld.wolfram.com/ContinuousFunction.html>`_.

A function is :math:`C^k` if the function, and every derivative of the function up to the :math:`k\text{th}` derivative is continuous.

We define the :math:`L^p` norm by: 

.. math::

   \lVert f \rVert_{L^p(\Omega)} := \Bigg( \int_\Omega | f | ^p \;dx \Bigg)^{\frac{1}{p}}

A function is :math:`L^p` integrable or in the :math:`L^p` space if :math:`\lVert f \rVert_{L^p(\Omega)} < \infty`

The Sobolev space :math:`V^{k, p}` is defined as 

.. math::

   V^{s, p}(\Omega) := \{ v \in L^p(\Omega) : \forall | \alpha | \leq s, \partial^\alpha_x v \in L^p(\Omega) \}

Here :math:`\alpha` is a multi index.
So this is saying for all of the partial derivatives of f who's total order are less than or equal to :math:`s`, that derivative is in :math:`L^p`.
When :math:`p` is left unspecified, we assume :math:`p=2`.

The space :math:`H^1` has :math:`L^2` integrable functions and first derivatives, and is defined as: 

.. math::

   H^1(\Omega) := V^{1, 2}(\Omega)

Using these spaces we define a few spaces to solve the vector-valued PDE with :math:`m` components on :math:`\Omega \subset \mathbb{R}^d`.

The first is used to represent the PDE solution by functions that are continuous internal to elements but can be discontinuous between elements. 
This **broken Soboloev space** is used for Discontinuous Galerkin methods.
Here, :math:`v\rvert_{\mathcal{K}}` is the function restricted to the domain of the element :math:`\mathcal{K}`.

.. math::

   V_u := \{ v \in [L^2(\Omega)]^m \;|\; \forall \mathcal{K} \in \mathcal{T},\; v\rvert_\mathcal{K} \in [H^1(\Omega)]^m \}

The next is another broken Sobolev space for primal formulations (explored in detail by Arnold et al [Arnold2000]_).

.. math::

   
   V_\sigma := \{ v \in [L^2(\Omega)]^{m\times d_x} \;|\; \forall \mathcal{K} \in \mathcal{T},\; v\rvert_\mathcal{K} \in [H^1(\Omega)]^{m \times d_x} \}

Where :math:`d_x = d-1` only considers the spatial dimensions.
Also define :math:`W_u` and :math:`W_\sigma` as the single valued trace spaces of :math:`V_u` and :math:`V_\sigma` respectively.

Next is a continuous function space over the domain, which will be used to represent the PDE solution, 
or continuous test function space for the interface conservation. 

.. math::

   V_{c} = \{ v \in [H^1(\Omega)]^m \}

and a corresponding space for primal formulations 

.. math::

   V_{p} = \{ v \in [H^1(\Omega)]^{m \times d_x} \}

And finally a continuous function space that will be used to aid in discretizing the mesh.

.. math:: 

   V_{y} = \{ v \in [H^1(\Omega)]^d \}

============================
Inner Products
============================

Weak formulations are composed out of :math:`L^2` inner products on the function spaces. 
We define the following notations for the inner products for scalars :math:`a, b` and vectors :math:`\mathbf{a}, \mathbf{b}`:

.. math::
   \langle a, b \rangle_\mathcal{K} &:= \sum_{\mathcal{K}\in\mathcal{T}} \int_\mathcal{K} a b \; d\mathcal{K} \\
   \langle \mathbf{a}, \mathbf{b} \rangle_\mathcal{K} &:= \sum_{\mathcal{K}\in\mathcal{T}} \int_\mathcal{K} \mathbf{a}\cdot\mathbf{b} \; d\mathcal{K} \\
   \langle a, b \rangle_\Gamma &:= \sum_{\Gamma\in\mathcal{E}} \int_\mathcal{\Gamma} a b \; d\Gamma \\
   \langle \mathbf{a}, \mathbf{b} \rangle_\Gamma &:= \sum_{\Gamma\in\mathcal{E}} \int_\Gamma \mathbf{a}\cdot\mathbf{b} \; d\Gamma


============================
MDG-ICE Formulation
============================

Consider the the conservation law:

.. math::

   \nabla \cdot \mathcal{F}(\mathbf{U}) = 0

Where the conservative flux :math:`\mathcal{F} : \mathbb{R}^m \mapsto \mathbb{R}^{m \times d}` can be split into convective and viscous components as follows:

.. math::

   \nabla \cdot \mathcal{F}^c(\mathbf{U}) + \nabla \cdot \mathcal{F}^v(\mathbf{U}, \nabla\mathbf{U}) = 0 \\
   \mathcal{F}^c = \begin{pmatrix}
      F^c_j \\ 
      \mathbf{U}
   \end{pmatrix}, 
   \quad \mathcal{F}^v = \begin{pmatrix} G_{ikrs}(U)\frac{\partial U_r}{\partial x_s} \\ 0 \end{pmatrix}

Here, the non-calligraphic :math:`F^c : \mathbb{R}^{m} \mapsto \mathbb{R}^{m \times d_x}` represents the spatial components of the convective flux. 
:math:`G_{ikrs}(\mathbf{U}) \in \mathbb{R}^{m\times d_x \times m \times d_x}` is the homogeneity tensor of the viscous flux, implying the viscous flux must be homogeneous with respect to the conservative variable gradients.
Also note that :math:`k` and :math:`s` are indices over the physical dimensions and omit the time dimension.

----------------------------
Original MDG-ICE Formulation
----------------------------
The original formulation recasts the conservation law into a system of first order equations with an auxiallary variable :math:`\sigma \in V_\sigma`.

.. math::

   \mathcal{F}(\mathbf{U}, \sigma) = \nabla\cdot \mathcal{F}^c + \nabla \cdot \begin{pmatrix} \sigma \\ 0 \end{pmatrix} &= 0 \text{ in } \mathcal{K} \quad&\forall \mathcal{K} \in \mathcal{T}\\
   \sigma - G_{ikrs}\frac{\partial U_r}{\partial x_s} &= 0 \text{ in } \mathcal{K} &\forall \mathcal{K} \in \mathcal{T} \\
   [[ n\cdot \mathcal{F}(\mathbf{U}, \sigma) ]] &= 0 \text{ on } \Gamma &\forall \Gamma \in \mathcal{E} \\
   \overline{G}_{ikrs}(U) [[ U_r n_s ]] &= 0 \text{ on } \Gamma &\forall \Gamma \in \mathcal{E} \\

where :math:`n` is the outward unit normal vector to the edge :math:`\Gamma`, :math:`\overline{u} = \frac{u_L + u_R}{2}` is the average operator,
and :math:`[[u]] = u_R - u_L` is the jump operator.


This results in the following weak form:

.. math::

   \langle \nabla \cdot \mathcal{F}(\mathbf{U}, \sigma), v_u \rangle_\mathcal{K}
   + \langle \sigma - G_{ikrs}\frac{\partial U_r}{\partial x_s}, v_\sigma \rangle_\mathcal{K} \\
   - \langle [[ n\cdot \mathcal{F}(\mathbf{U}, \sigma) ]], w_u \rangle_\Gamma
   - \langle \overline{G}_{ikrs}(U) [[ U_r n_s ]], w_\sigma \rangle_\Gamma \\
   = 0  \quad \forall v_u\in V_u, v_\sigma \in V_\sigma, w_u \in W_u, w_\sigma \in W_\sigma

The domain integrals :math:`\langle \cdot, \cdot \rangle_\mathcal{K}` can be integrated by parts and resulting fluxes in the trace integrals replaced by consistent and conservative numerical fluxes.

-------------------- 
Proposed Formulation
--------------------

.. math:: 
   \mathcal{F}(\mathbf{U}, \nabla \mathbf{U}) &= 0 \text{ in } \mathcal{K} \quad&\forall \mathcal{K} \in \mathcal{T}\\
   [[ n\cdot \mathcal{F}(\mathbf{U}, \sigma) ]] &= 0 \text{ on } \Gamma &\forall \Gamma \in \mathcal{E} \\
   \overline{G}_{ikrs}(U) [[ U_r n_s ]] &= 0 \text{ on } \Gamma &\forall \Gamma \in \mathcal{E} \\

We want to also investigate if the last set of equations (the interface conservation for the constitutive relation) is necessary. 
The corresponding weak form is:

.. math:: 

   \langle \nabla \cdot \mathcal{F}(\mathbf{U}, \nabla \mathbf{U}), v_u \rangle_\mathcal{K} \\
   - \langle [[ n\cdot \mathcal{F}(\mathbf{U}, \sigma) ]], v_c \rangle_\Gamma
   - \langle \overline{G}_{ikrs}(U) [[ U_r n_s ]], v_p \rangle_\Gamma \\
   = 0  \quad \forall v_u\in V_u, v_c \in V_c, v_p \in V_p

-------------------------
Geometry Parameterization
-------------------------

For the purpose of the discretization the geometry (described by node coordinates) exists in the :math:`V_y` space. However, this neglects the enforcement of boundary conditions.
Corrigan et al. introduce an operator :math:`b : V_y \mapsto V_y` that enforces the boundary condition.
We, instead parameterize the geometry by a function :math:`g : \mathbb{R}^{n_p} \mapsto V_y`. 
Where :math:`n_p` is the number of parameters for the parameterization. 
This also allows for parameterizations like setting the final time by a single input variable.

----------------------
Variational Derivative
----------------------

Define the variational derivative of a functional :math:`\Pi(u)` as: 

.. math::

   \Pi_u(u; \hat{u}) := \frac{d}{d\epsilon} \Pi(u + \epsilon \hat{u}) \Bigr|_{\epsilon = 0}

---------------------
Regularization
---------------------

Let :math:`\tilde{u} \in \mathbb{R}^{n_p}` represent the geometry parameterization.
Isotropic Laplacian-type mesh regularization results in the following operator:

.. math::

   -\langle \nabla g_\tilde{u}(\tilde{u}; \hat{u}), \nabla g_\tilde{u}(\tilde{u}; \hat{v}) \rangle_\mathcal{K}

If we denote the element stiffness matrix by 

.. math::

   K_\mathcal{K} = \int_\mathcal{K} \nabla u \cdot \nabla v \; d\mathcal{K}

Then this isotropic Laplacian-type mesh regularization can be written as:

.. math::

   -(g_\tilde{u}(\tilde{u}; \hat{u})\mathbf{I})^T K_\mathcal{K} g_\tilde{u}(\tilde{u}; \hat{v}) \mathbf{I}

which can be implemented as a simple loop over the elements.

==========
References
==========
.. [Kercher2021] Kercher, A. D., Corrigan, A., & Kessler, D. A. (2021). The moving discontinuous Galerkin finite element method with interface condition enforcement for compressible viscous flows. International Journal for Numerical Methods in Fluids, 93(5), 1490-1519.

.. [Arnold2000]  Arnold, D. N., Brezzi, F., Cockburn, B., & Marini, D. (2000). Discontinuous Galerkin methods for elliptic problems. In Discontinuous Galerkin Methods: Theory, Computation and Applications (pp. 89-101). Berlin, Heidelberg: Springer Berlin Heidelberg.

.. [GrundmannMoller1978] Grundmann, A., & Möller, H. M. (1978). Invariant integration formulas for the n-simplex by combinatorial methods. SIAM Journal on Numerical Analysis, 15(2), 282-290.
