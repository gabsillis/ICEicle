local gamma = 1.4
local alpha = 1
local K = 5.0
local rho_inf = 1.0
local p_inf = 1.0 / gamma
local u_inf = 2
local v_inf = 2
local x0 = -10
local y0 = -10
local Rgas = 287.052874

local umag2_inf = u_inf * u_inf + v_inf * v_inf
local rhoE_inf = p_inf / (gamma - 1) + 0.5 * rho_inf * umag2_inf
local csound_inf = math.sqrt(gamma * p_inf / rho_inf)
local T_inf = p_inf / (rho_inf * Rgas)

local tfinal = 5

local exact = function(t)
	return function(x, y)
		local xbar = x - x0 - u_inf * t;
		local ybar = y - y0 - v_inf * t;
		local rbar = math.sqrt(xbar ^ 2 + ybar ^ 2)

		local u = u_inf - K / (2 * math.pi) * ybar * math.exp(alpha * (1 - rbar ^ 2) / 2)
		local v = v_inf + K / (2 * math.pi) * xbar * math.exp(alpha * (1 - rbar ^ 2) / 2)
		-- T / T_inf
		local T_hat = 1 -
			K ^ 2 * (gamma - 1) / (8 * alpha * math.pi ^ 2 * csound_inf ^ 2) * math.exp(alpha * (1 - rbar ^ 2))

		local rho = rho_inf * (T_hat) ^ (1 / (gamma - 1))
		local p = rho_inf * (T_hat) ^ (gamma / (gamma - 1))

		local umag2 = u * u + v * v
		local rhoE = p / (gamma - 1) + 0.5 * rho * umag2
		return { rho, rho * u, rho * v, rhoE }
	end
end

return {
	ndim = 2,
	uniform_mesh = {
		nelem = { 20, 20 },
		bounding_box = { min = { -20, -20 }, max = { 20, 20 } },
		boundary_conditions = {
			types = {
				"riemann",
				"extrapolation",
				"extrapolation",
				"extrapolation"
			},
			flags = { 0, 0, 0, 0 },
			geometry_order = 1,
		},
	},

	-- define the finite element domain
	fespace = {
		-- the basis function type (optional: default = lagrange)
		basis = "lagrange",

		-- the quadrature type (optional: default = gauss)
		quadrature = "gauss",

		-- the basis function order
		order = 1,

	},

	initial_condition = exact(0),

	conservation_law = {
		name = "euler",
		free_stream = { rho_inf, rho_inf * u_inf, rho_inf * v_inf, rhoE_inf },
	},

	solver = {
		type = "rk3-ssp",
		cfl = 0.1,
		tfinal = tfinal,
		ivis = 1,
		idiag = 1,
	},

	output = {
		writer = "vtu"
	},

	post = {
		exact_solution = exact(tfinal),
		tasks = {
			"l2_error",
			"l1_error",
			"plot_exact_projection"
		},
	}


}
