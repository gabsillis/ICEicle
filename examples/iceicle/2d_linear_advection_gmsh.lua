-- Example: Time Dependent Linear Advection in 2D
-- Author Gianni Absillis
--
-- Solves a normal distribution bump being transported by linear advection
-- with a triangular mesh read in from gmsh

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

    -- read in a mesh from gmsh
    gmsh = {
        file = "../gmsh_meshes/square_2x2.msh",

        bc_definitions = {
            -- 1: airfoil boundary
            { "dirichlet", 0 },

            -- 2: bottom and top walls
            { "dirichlet", 0 },

            -- 3: outlet
            { "dirichlet", 0 },

            -- 4: inlet
            { "dirichlet", 0 },
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
