-- This is a test problem from "I do like CFD"
-- This 1D problem forms a viscous shock structure

local gamma = 1.4

-- left state
local rhoL = 1.0
local machL = 0.6
local pL = 1 / gamma
local TL = gamma * pL / rhoL

-- right state
local rhoR = (gamma + 1) * machL ^ 2 / ((gamma - 1) * machL ^ 2 + 2)
local uR = rhoL * machL / rhoR
local pR = 1 / gamma * (1 + (2 * gamma) / (gamma + 1) * (machL ^ 2 - 1))
local TR = gamma * pR / rhoR
