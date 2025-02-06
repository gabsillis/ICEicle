return {
    ndim = 2,

    uniform_mesh = {
        nelem = { 32, 32 },
        bounding_box = {
            min = { 0.0, 0.0 },
            max = { 1.0, 1.0 },
        },

        -- quad_ratio = { 0.0, 0.0 },

        -- set boundary conditions
        boundary_conditions = {
            -- the boundary condition types
            -- in order of direction and side
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

    -- set up Navier Stokes
    conservation_law = {
        name = "navier-stokes",

        viscosity = {
            name = "constant",
            mu = 0.1,
        },

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
    },

    -- boundary conditions
    boundary_conditions = {
        dirichlet = {
            function(x, y)
                local rho = math.sin(2 * (x + y)) + 4
                -- local rho = 1 + x
                return rho, rho * 1, rho * 1, rho ^ 2
            end,
        },
    },

    -- initial condition
    initial_condition = function(x, y)
        return {
            math.sin(1.5 * (x + y)) + 3,
            0.2 * math.sin(1.5 * (x + y)) + 3,
            0.8 * math.sin(1.5 * (x + y)) + 3,
            (math.sin(1.5 * (x + y)) + 3) ^ 2
        }
    end,

    -- solver
    solver = {
        type = "gauss-newton",
        form_subproblem_mat = true,
        lambda_u = 1e-8,
        kmax = 100,
        tau_abs = 1e-10,
        ivis = 1,
        -- idiag = 1,
    },

    -- output
    output = {
        writer = "vtk",
    },

    -- post-processing
    post = {
        exact_solution = function(x, y)
            local rho = math.sin(2 * (x + y)) + 4
            -- local rho = 1 + x
            return { rho, rho * 1, rho * 1, rho ^ 2 }
        end,

        tasks = {
            "l2_error",
            "plot_exact_projection"
        },
    },
}
