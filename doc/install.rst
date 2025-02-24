Installation
============

==========
Quickstart 
==========
A bare-bones version of the library (with unit tests built by default) can be built with 

.. code-block:: bash

   mkdir build 
   cd build 
   cmake ..
   make -j

=====================
Third Party Libraries
=====================

Third party libraries can be included in the build process by applying the corresponding cmake flag:
``cmake .. -D<flag-name>=ON``. 
Unless otherwise stated, the library is automatically downloaded using FetchContent by the cmake configuration when the flag is enabled.

---
MPI   
---
``ICEICLE_USE_MPI``: is on by default and compiles and links with the system MPI for parallel domain decomposition

-----
PETSc
-----
``ICEICLE_USE_PETSC``: Interface ICEicle with PETSc - a library for linear and nonlinear sparse solvers on distributed
memory platforms. Currently integrated features include:
- Finite difference construction of jacobians in a PETSc ``Mat``.
- Newton solver that uses ``KSP`` linear solvers 
- Gauss-Newton solver that uses ``KSP`` linear solvers 
- Matrix free Newton-Krylov solver 

This library is located by use of the ``PETSC_DIR`` and ``PETSC_ARCH`` environment variables using pkgconfig.

---
Lua
---
``ICEICLE_USE_LUA``: Interface with thelua C api through sol2. This is used for input decks with callback functions.

.. note::
   This requires Lua development libraries.

To get the lua development libraries on Fedora or RHEL

.. code-block:: bash

   dnf install lua-devel

-----
METIS
-----
``ICEICLE_USE_METIS``: Enable mesh partitioning with METIS. The petsc bitbucket fork is used for the automatic install.

---
VTK
---

``ICEICLE_USE_VTK`` interface with the VTK c++ api. This is needed for parallel output in greater than 1 dimension.

.. note::
   This requires a VTK built with mpi. It may be necessary to install VTK from the source and specify ``VTK_DIR``

-------------
Documentation 
-------------

These documentation pages can be built by enabling ``ICEICLE_BUILD_DOC``. This requires Sphinx and doxygen to be installed, along with the packages ``breathe``, ``sphinx_rtd_theme``, ``sphinx.ext.mathjax``, ``sphinx.ext.autosectionlabel``, and ``sphinxcontrib.tikz``.
