-- Example: Gmsh Spacetime Burgers equation
-- Author: Gianni Absillis (gabsill@ncsu.edu)
--
-- This example solves the burgers equation on a spacetime domain given by
-- a gmsh mesh generated by st_burgers.py in gmsh_meshes

local t_s = 0.5

local y_inf = 0.2

-- return the configuration table
return {
	-- specify the number of dimensions (REQUIRED)
	ndim = 2,

	-- read in a mesh from gmsh
	gmsh = {
		file = "../gmsh_meshes/st_burgers.msh",
		bc_definitions = {
			-- 1: initial condition (dirichlet, flag 1)
			{ "dirichlet",        1 },

			-- 2: right (extrapolation, flag 0)
			{ "extrapolation",    0 },

			-- 3: spacetime-future
			{ "spacetime-future", 0 },

			-- 4: left (dirichlet, flag 0)
			{ "dirichlet",        0 },
		},
	},

	-- define the finite element domain
	fespace = {
		-- the basis function type (optional: default = lagrange)
		basis = "lagrange",

		-- the quadrature type (optional: default = gauss)
		quadrature = "gauss",

		-- the basis function order
		order = 2,
	},

	-- describe the conservation law
	conservation_law = {
		-- the name of the conservation law being solved
		name = "spacetime-burgers",
		mu = 1e-3,
		a_adv = { 0.0 },
		b_adv = { 1.0 },
	},

	-- initial condition
	initial_condition = function(x, t)
		-- return y_inf
		return 1.0 / (2 * math.pi * t_s) * math.sin(2 * math.pi * x) + y_inf
	end,

	-- boundary conditions
	boundary_conditions = {
		dirichlet = {
			y_inf,
			function(x, t)
				-- return y_inf
				return 1.0 / (2 * math.pi * t_s) * math.sin(2 * math.pi * x) + y_inf
				-- return math.exp(-0.5 * ((x - 0.5) / 0.1) ^ 2) / 0.1;
			end,
		},
	},

	-- MDG
	--	mdg = {
	--		ncycles = 1,
	--		ic_selection_threshold = function(icycle)
	--			return 0
	--		end,
	--	},

	-- solver
	solver = {
		type = "mfnk",
		form_subproblem_mat = false,
		linesearch = {
			type = "none",
		},
		lambda_b = 1e-1,
		lambda_u = 0,
		ivis = 1,
		tau_abs = 1e-10,
		tau_rel = 0,
		kmax = 10,
	},

	-- output
	output = {
		writer = "vtu",
	},
}
