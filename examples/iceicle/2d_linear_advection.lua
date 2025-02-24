-- Example: Time Dependent Linear Advection in 2D
-- Author Gianni Absillis
--
-- Solves a normal distribution bump being transported by linear advection

-- PDF of the normal distribution centered at (mu_x, mu_y) with standandard deviations sigmax, sigmay
local normal_dist = function(mu_x, mu_y, sigma_x, sigma_y)
    return function(x, y)
        return 1.0 / (2 * math.pi * sigma_x * sigma_y) * math.exp(
            -1.0 / 2.0 * (
                ((x - mu_x) / sigma_x) ^ 2
                + ((y - mu_y) / sigma_y) ^ 2
            )
        );
    end
end

-- This is the table that gets read in by the program
-- and will be processed
return {
    -- sepecify the number of dimensions (REQUIRED)
    ndim = 2,

    -- create a uniform mesh
    uniform_mesh = {
        -- specify the number of elements in each direction
        -- (there is only the x direction in 1D)
        nelem = { 100, 100 },

        -- Specify a bounding hypercube in R^ndim space
        -- using the extreme corners
        bounding_box = {
            min = { 0.0, 0.0 },
            max = { 2.0, 2.0 },
        },

        -- set the boundary conditions of the sides of the bounding box
        boundary_conditions = {
            -- the boundary condition types
            -- in order first by direction, then positive/negative side
            types = {
                "dirichlet", -- left side
                "dirichlet", -- bottom side
                "dirichlet", -- right side
                "dirichlet", -- top side
            },

            -- the boundary condition flags
            -- used to identify user defined state
            flags = {
                0, -- left
                0, -- bottom
                0, -- right
                0, -- top
            },
        },
        geometry_order = 1,
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
        -- in this case, a special case of burgers equation
        name = "burgers",

        -- the linear advection speed
        a_adv = { 1.0, 0.5 },

        -- the diffusion coefficient
        mu = 0.0,
    },

    -- initial condition: create a gauss hump centered at (0.5, 0.5)
    initial_condition = normal_dist(0.5, 0.5, 0.1, 0.1),

    -- specify boundary conditions
    boundary_conditions = {
        dirichlet = {
            0.0, -- flag 0
        }
    },

    -- solve with an explicit solver
    solver = {
        type = "rk3-tvd", -- the type of the solver
        cfl = 0.1,
        -- ntime = 1,
        tfinal = 1.0,
        ivis = 100,
    },

    -- output
    output = {
        writer = "vtk"
    },

    -- post-processing
    post = {
        -- exact solution is the initial condition moved over at t=1
        exact_solution = normal_dist(1.5, 1.0, 0.1, 0.1),

        -- we want to report the L2 function norm of the error
        -- so we add this to the tasks
        tasks = {
            "l1_error",
            "l2_error",
        },
    },
}
