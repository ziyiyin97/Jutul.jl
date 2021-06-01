export simulate, perform_step!
export Simulator, TervSimulator
export simulator_config

abstract type TervSimulator end
struct Simulator <: TervSimulator
    model::TervModel
    storage::NamedTuple
end

function Simulator(model; state0 = nothing, parameters = setup_parameters(model), copy_state = true, kwarg...)
    # We need to sort the secondary variables according to their dependency ordering before simulating.
    sort_secondary_variables!(model)
    if isnothing(state0)
        state0 = setup_state(model)
    elseif copy_state
        # Take a deep copy to avoid side effects.
        state0 = deepcopy(state0)
    end
    storage = setup_storage(model, state0 = state0, parameters = parameters)
    # Initialize for first time usage
    initialize_storage!(storage, model; kwarg...)
    # We convert the mutable storage (currently Dict) to immutable (NamedTuple)
    # This allows for much faster lookup in the simulation itself.
    storage = convert_to_immutable_storage(storage)
    Simulator(model, storage)
end

function perform_step!(simulator::TervSimulator; vararg...)
    perform_step!(simulator.storage, simulator.model; vararg...)
end

function perform_step!(storage, model; dt = nothing, linsolve = nothing, forces = nothing, iteration = NaN, config = simulator_config(sim))
    timing_out = config[:debug_level] > 1
    # Update the properties and equations
    t_asm = @elapsed begin 
        update_state_dependents!(storage, model, dt, forces)
    end
    if timing_out
        @debug "Assembled equations in $t_asm seconds."
    end
    # Update the linearized system
    t_lsys = @elapsed begin
        update_linearized_system!(storage, model)
    end
    if timing_out
        @debug "Updated linear system in $t_lsys seconds."
    end
    converged, e, tol = check_convergence(storage, model, iteration = iteration, dt = dt, extra_out = true)
    if converged
        do_solve = iteration == 1
        @debug "Step converged."
    else
        do_solve = true
    end
    if do_solve
        t_solve, t_update = solve_update!(storage, model::TervModel; linsolve = linsolve)
        if timing_out
            @debug "Solved linear system in $t_solve seconds."
            @debug "Updated state $t_update seconds."
        end
    end
    return (e, tol)
end

function simulator_config(sim)
    cfg = Dict()
    cfg[:max_timestep_cuts] = 5
    cfg[:max_nonlinear_iterations] = 10
    cfg[:linear_solver] = nothing
    cfg[:output_states] = true
    # Define debug level. If debugging is on, this determines the amount of output.
    cfg[:debug_level] = 1
    return cfg
end

function simulate(sim::TervSimulator, timesteps::AbstractVector; forces = nothing, config = simulator_config(sim))
    states = []
    no_steps = length(timesteps)
    maxIterations = config[:max_nonlinear_iterations]
    linsolve = config[:linear_solver]
    @info "Starting simulation"
    for (step_no, dT) in enumerate(timesteps)
        t_str =  Dates.canonicalize(Dates.CompoundPeriod(Second(dT)))
        @info "Solving step $step_no/$no_steps of length $t_str."
        dt = dT
        done = false
        t_local = 0
        cut_count = 0
        while !done
            ok = solve_ministep(sim, dt, maxIterations, linsolve, forces, config)
            if ok
                t_local += dt
                if t_local >= dT
                    break
                end
            else
                max_cuts = config[:max_timestep_cuts]
                if cut_count > max_cuts
                    @warn "Unable to converge time step $step_no/$no_steps. Aborting."
                    return states
                end
                dt = min(dt/2, dT - t_local)
                @warn "Cutting time-step. Step $(100*t_local/dT) % complete.\nStep fraction reduced to $(100*dt/dT)% of full step.\nThis is cut $cut_count of $max_cuts allowed."
                cut_count += 1
            end
        end
        if config[:output_states]
            store_output!(states, sim)
        end
    end
    return states
    @info "Simulation complete."
end

function solve_ministep(sim, dt, maxIterations, linsolve, forces, cfg)
    done = false
    for it = 1:maxIterations
        e, tol = perform_step!(sim, dt = dt, iteration = it, forces = forces, linsolve = linsolve, config = cfg)
        done = e < tol
        if done
            break
        end
        if e > 1e10 || isinf(e) || isnan(e)
            break
        end
    end

    if done
        t_finalize = @elapsed update_after_step!(sim)
        if cfg[:debug_level] > 1
            @debug "Finalized in $t_finalize seconds."
        end
    else
        primary = sim.storage.primary_variables
        for f in keys(primary)
            update_values!(primary[f], sim.storage.state0[f])
        end
    end
    return done
end

function update_after_step!(sim)
    update_after_step!(sim.storage, sim.model)
end

function store_output!(states, sim)
    state_out = get_output_state(sim.storage, sim.model)
    push!(states, state_out)
end
