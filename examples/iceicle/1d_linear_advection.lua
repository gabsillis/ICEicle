-- Example: Time Dependent Linear Advection
-- Author: Gianni Absillis
--
-- Solves a gauss hump being transported by linear advection

-- PDF of the normal distribution, forms a smooth bump
-- centered at mean with standard deviation stddev
local normal_dist = function(mean, stddev)
    return function(x)
        return 1.0 / math.sqrt(2 * math.pi * stddev ^ 2) * math.exp(-(x - mean) ^ 2 / (2 * stddev ^ 2));
    end
end

-- This is the table that gets read in by the program
-- and will be processed
return {

    -- sepecify the number of dimensions (REQUIRED)
    ndim = 1,

    -- create a uniform mesh
    uniform_mesh = {

        -- specify the number of elements in each direction
        -- (there is only the x direction in 1D)
        nelem = { 100 },

        -- Specify a bounding hypercube in R^ndim space
        -- using the extreme corners
        bounding_box = {
            min = { 0.0 },
            max = { 2.0 },
        },

        -- set the boundary conditions of the sides of the bounding box
        boundary_conditions = {

            -- the boundary condition types
            -- in order first by direction, then positive/negative side
            types = {
                "dirichlet", -- left side
                "dirichlet", -- right side
            },

            -- the boundary condition flags
            -- used to identify user defined state
            flags = {
                0, -- left
                1, -- right
            },
        }
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
        a_adv = { 1.0 },

        -- the diffusion coefficient
        mu = 0.0,
    },

    -- initial condition: create a gauss hump centered at 0.5
    initial_condition = normal_dist(0.5, 0.1),

    -- specify boundary conditions
    boundary_conditions = {
        dirichlet = {
            0.0, -- flag 0
            0.0, -- flag 1
        }
    },

    -- solve with an explicit solver
    solver = {
        type = "rk3-tvd", -- the type of the solver
        cfl = 0.1,
        tfinal = 1.0,
        ivis = 1000,
    },

    -- output
    output = {
        writer = "dat",
    },

    -- post-processing
    post = {

        -- exact solution is the initial condition moved over by 1.0
        exact_solution = normal_dist(1.5, 0.1),

        -- we want to report the L2 function norm of the error
        -- so we add this to the tasks
        tasks = {
            "l1_error",
            "plot_exact_projection", -- we want to plot a projection of the exact solution
        },
    },

}
