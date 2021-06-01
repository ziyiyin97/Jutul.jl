export simulate, perform_step!
export Simulator, TervSimulator


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

function perform_step!(storage, model; dt = nothing, linsolve = nothing, forces = nothing, iteration = NaN)
    # Update the properties and equations
    update_state_dependents!(storage, model, dt, forces)
    # Update the linearized system
    t_lsys = @elapsed begin
        update_linearized_system!(storage, model)
    end
    @debug "Updated linear system in $t_lsys seconds."

    converged, e, tol = check_convergence(storage, model, iteration = iteration, dt = dt, extra_out = true)
    if converged
        do_solve = iteration == 1
        @debug "Step converged."
    else
        do_solve = true
    end
    if do_solve
        t_solve, t_update = solve_update!(storage, model::TervModel; linsolve = linsolve)
        @debug "Solved linear system in $t_solve seconds."
        @debug "Updated state $t_update seconds."
    end
    return (e, tol)
end

function simulate(sim::TervSimulator, timesteps::AbstractVector; maxIterations = 10, outputStates = true, forces = nothing, linsolve = nothing)
    states = []
    no_steps = length(timesteps)
    @info "Starting simulation"
    for (step_no, dT) in enumerate(timesteps)
        t_str =  Dates.canonicalize(Dates.CompoundPeriod(Second(dT)))
        @info "Solving step $step_no/$no_steps of length $t_str."
        dt = dT
        done = false
        t_local = 0
        cut_count = 0
        while !done
            ok = solve_ministep(sim, dt, maxIterations, linsolve, forces)
            if ok
                t_local += dt
                if t_local >= dT
                    break
                end
            else
                @warn "Cutting time-step."
                @assert cut_count < 5
                dt = min(dt/2, dT - t_local)
                cut_count += 1
            end
        end
        if outputStates
            store_output!(states, sim)
        end
    end
    return states
    @info "Simulation complete."
end

function solve_ministep(sim, dt, maxIterations, linsolve, forces)
    done = false
    for it = 1:maxIterations
        e, tol = perform_step!(sim, dt = dt, iteration = it, forces = forces, linsolve = linsolve)
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
        @debug "Finalized in $t_finalize seconds."
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
