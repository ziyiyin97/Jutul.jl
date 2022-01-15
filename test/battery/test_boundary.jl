#=
Simple current collector
A conductro with constant conductivity
=#
using Revise
using Jutul
using Test
#include("../../src/applications/battery/models/CurrentAndVoltageBoundary.jl")
#ENV["JULIA_DEBUG"] = Jutul;



  
timesteps = [10.,10.]


# System type for function overloading
sys = CurrentAndVoltageSystem()
domain = CurrentAndVoltageDomain()
    # Setup model
model = SimulationModel(domain, sys, context = DefaultContext())

# State is dict with pressure in each cell
phi = 1.0
current = 1.0
boudary_phi = [1., 2.]
#S = model.secondary_variables
#S[:Phi] = Volt()
#S[:Current] = Current()

# Inital values for variable. Variables w/o update_as_secondary must be set here
init = Dict(:Phi => phi, :Current=>current)
state0 = setup_state(model, init)
        
# Model parameters
parameters = setup_parameters(model)
parameters[:tolerances][:default] = 1e-8

# Contains storage
sim = Simulator(model, state0=state0, parameters=parameters)

cfg = simulator_config(sim)
cfg[:linear_solver] = nothing

# Run simulation
currentFun(time) = time
forces = Dict( :current => currentFun)
states, _ = simulate(sim, timesteps, config = cfg, forces = forces)

##

