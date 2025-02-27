export state_gradient, solve_adjoint_sensitivities, solve_adjoint_sensitivities!, setup_adjoint_storage

"""
    solve_adjoint_sensitivities(model, states, reports, G; extra_timing = false, state0 = setup_state(model), forces = setup_forces(model), raw_output = false, kwarg...)

Compute sensitivities of `model` parameter with name `target` for objective function `G`.

Solves the adjoint equations: For model equations F the gradient with respect to parameters p is
    ∇ₚG = Σₙ (∂Fₙ / ∂p)ᵀ λₙ where n ∈ [1, N].
Given Lagrange multipliers λₙ from the adjoint equations
    (∂Fₙ / ∂xₙ)ᵀ λₙ = - (∂J / ∂xₙ)ᵀ - (∂Fₙ₊₁ / ∂xₙ)ᵀ λₙ₊₁
where the last term is omitted for step n = N and G is the objective function.
"""
function solve_adjoint_sensitivities(model, states, reports, G; 
                                        n_objective = nothing,
                                        extra_timing = false,
                                        state0 = setup_state(model),
                                        forces = setup_forces(model),
                                        raw_output = false,
                                        extra_output = false, kwarg...)
    # One simulator object for the equations with respect to primary (at previous time-step)
    # One simulator object for the equations with respect to parameters
    set_global_timer!(extra_timing)
    # Allocation part
    storage = setup_adjoint_storage(model; state0 = state0, n_objective = n_objective, kwarg...)
    parameter_model = storage.parameter.model
    n_param = number_of_degrees_of_freedom(parameter_model)
    ∇G = gradient_vec_or_mat(n_param, n_objective)
    # Timesteps
    timesteps = report_timesteps(reports)
    N = length(states)
    @assert length(reports) == N == length(timesteps)
    # Solve!
    solve_adjoint_sensitivities!(∇G, storage, states, state0, timesteps, G, forces = forces)
    print_global_timer(extra_timing; text = "Adjoint solve detailed timing")
    if raw_output
        out = ∇G
    else
        out = store_sensitivities(parameter_model, ∇G, storage.parameter_map)
    end
    if extra_output
        return (out, storage)
    else
        return out
    end
end

"""
    setup_adjoint_storage(model; state0 = setup_state(model), parameters = setup_parameters(model))

Set up storage for use with `solve_adjoint_sensitivities!`.
"""
function setup_adjoint_storage(model; state0 = setup_state(model),
                                      parameters = setup_parameters(model),
                                      n_objective = nothing,
                                      targets = parameter_targets(model),
                                      use_sparsity = true,
                                      linear_solver = select_linear_solver(model, mode = :adjoint, rtol = 1e-8),
                                      param_obj = false, kwarg...)
    primary_model = adjoint_model_copy(model)
    # Standard model for: ∂Fₙᵀ / ∂xₙ
    forward_sim = Simulator(primary_model, state0 = deepcopy(state0), parameters = deepcopy(parameters), mode = :forward, extra_timing = nothing)
    # Same model, but adjoint for: ∂Fₙ₊₁ᵀ / ∂xₙ
    backward_sim = Simulator(primary_model, state0 = deepcopy(state0), parameters = deepcopy(parameters), mode = :reverse, extra_timing = nothing)
    # Create parameter model for ∂Fₙ / ∂p
    parameter_model = adjoint_parameter_model(model, targets)
    # Note that primary is here because the target parameters are now the primaries for the parameter_model
    parameter_map, = variable_mapper(parameter_model, :primary, targets = targets; kwarg...)
    # Transfer over parameters and state0 variables since many parameters are now variables
    state0_p = swap_variables(state0, parameters, parameter_model, variables = true)
    parameters_p = swap_variables(state0, parameters, parameter_model, variables = false)
    parameter_sim = Simulator(parameter_model, state0 = deepcopy(state0_p), parameters = deepcopy(parameters_p), mode = :sensitivities, extra_timing = nothing)
    if use_sparsity isa Bool
        if use_sparsity
            # We will update these later on
            sparsity_obj = Dict{Any, Any}(:parameter => nothing, 
                                        :forward => nothing)
        else
            sparsity_obj = nothing
        end
    else
        # Assume it was manually set up
        sparsity_obj = use_sparsity
    end
    n_pvar = number_of_degrees_of_freedom(model)
    λ = gradient_vec_or_mat(n_pvar, n_objective)
    fsim_s = forward_sim.storage
    rhs = vector_residual(fsim_s.LinearizedSystem)
    dx = fsim_s.LinearizedSystem.dx_buffer
    n_var = length(dx)
    if param_obj
        n_param = number_of_degrees_of_freedom(parameter_model)
        dobj_dparam = gradient_vec_or_mat(n_param, n_objective)
        param_buf = similar(dobj_dparam)
    else
        dobj_dparam = nothing
        param_buf = nothing
    end
    rhs_transfer_needed = length(rhs) != n_var
    multiple_rhs = !isnothing(n_objective)
    if multiple_rhs
        # Need bigger buffers for multiple rhs
        rhs = gradient_vec_or_mat(n_var, n_objective)
        dx = gradient_vec_or_mat(n_var, n_objective)
    elseif rhs_transfer_needed
        rhs = zeros(n_var)
    end
    return (forward = forward_sim,
            backward = backward_sim,
            parameter = parameter_sim,
            parameter_map = parameter_map,
            objective_sparsity = sparsity_obj,
            lagrange = λ,
            lagrange_buffer = similar(λ),
            dparam = dobj_dparam,
            param_buf = param_buf,
            dx = dx,
            rhs = rhs,
            linear_solver = linear_solver,
            multiple_rhs = multiple_rhs,
            rhs_transfer_needed = rhs_transfer_needed,
            n = number_of_degrees_of_freedom(parameter_model)
            )
end

"""
    solve_adjoint_sensitivities!(∇G, storage, states, state0, timesteps, G; forces = setup_forces(model))

Non-allocating version of `solve_adjoint_sensitivities`.
"""
function solve_adjoint_sensitivities!(∇G, storage, states, state0, timesteps, G; forces = setup_forces(model))
    N = length(timesteps)
    @assert N == length(states)
    # Do sparsity detection if not already done.
    update_objective_sparsity!(storage, G, states, timesteps, forces, :forward)
    update_objective_sparsity!(storage, G, states, timesteps, forces, :parameter)
    # Set gradient to zero before solve starts
    @. ∇G = 0
    @tic "sensitivities" for i in N:-1:1
        fn = deepcopy
        if i == 1
            s0 = fn(state0)
        else
            s0 = fn(states[i-1])
        end
        if i == N
            s_next = nothing
        else
            s_next = fn(states[i+1])
        end
        s = fn(states[i])
        update_sensitivities!(∇G, i, G, storage, s0, s, s_next, timesteps, forces)
    end
    dparam = storage.dparam
    if !isnothing(dparam)
        @. ∇G += dparam
    end
    rescale_sensitivities!(∇G, storage.parameter.model, storage.parameter_map)
    @assert all(isfinite, ∇G)
    return ∇G
end

function update_objective_sparsity!(storage, G, states, timesteps, forces, k = :forward)
    obj_sparsity = storage.objective_sparsity
    if isnothing(obj_sparsity) || (k == :parameter && isnothing(storage.dparam))
        return
    else
        sparsity = obj_sparsity[k]
        if isnothing(sparsity)
            sim = storage[k]
            obj_sparsity[k] = determine_objective_sparsity(sim, sim.model, G, states, timesteps, forces)
        end
    end
end

function get_objective_sparsity(storage, k)
    obj_sparsity = storage.objective_sparsity
    if isnothing(obj_sparsity)
        S = nothing
    else
        S = obj_sparsity[k]
    end
    return S
end

function determine_objective_sparsity(sim, model, G, states, timesteps, forces)
    update_secondary_variables!(sim.storage, sim.model)
    state = sim.storage.state
    # m, state, dt, step_no, forces
    F_outer = (state, i) -> G(model, state, timesteps[i], i, forces_for_timestep(sim, forces, timesteps, i))
    sparsity = determine_sparsity_simple(s -> F_outer(s, 1), model, state)
    for i in 2:length(states)
        s_new = determine_sparsity_simple(s -> F_outer(s, i), model, state)
        merge_sparsity!(sparsity, s_new)
    end
    return sparsity
end

function merge_sparsity!(sparsity::AbstractDict, s_new::AbstractDict)
    for (k, v) in s_new
        merge_sparsity!(sparsity[k], v)
    end
end
function merge_sparsity!(sparsity::AbstractVector, s_new::AbstractVector)
    for k in s_new
        push!(sparsity, k)
    end
    unique!(sparsity)
end

function state_gradient(model, state, F, extra_arg...; kwarg...)
    n_total = number_of_degrees_of_freedom(model)
    ∂F∂x = zeros(n_total)
    return state_gradient!(∂F∂x, model, state, F, extra_arg...; kwarg...)
end

function state_gradient!(∂F∂x, model, state, F, extra_arg...; parameters = setup_parameters(model))
    # Either with respect to all primary variables, or all parameters.
    state = setup_state(model, state)
    state = merge_state_with_parameters(model, state, parameters)
    state = convert_state_ad(model, state)
    state = convert_to_immutable_storage(state)
    update_secondary_variables_state!(state, model)
    state_gradient_outer!(∂F∂x, F, model, state, extra_arg)
    return ∂F∂x
end

function merge_state_with_parameters(model, state, parameters)
    for (k, v) in pairs(parameters)
        state[k] = v
    end
    return state
end

function state_gradient_outer!(∂F∂x, F, model, state, extra_arg; sparsity = nothing)
    if !isnothing(sparsity)
        @. ∂F∂x = 0
    end
    state_gradient_inner!(∂F∂x, F, model, state, nothing, extra_arg; sparsity = sparsity)
end

function state_gradient_inner!(∂F∂x, F, model, state, tag, extra_arg, eval_model = model; sparsity = nothing)
    layout = matrix_layout(model.context)
    get_partial(x::AbstractFloat, i) = 0.0
    get_partial(x::ForwardDiff.Dual, i) = x.partials[i]

    function store_partials!(∂F∂x::AbstractVector, v, i, ne, np, offset)
        for p_i in 1:np
            ix = alignment_linear_index(i, p_i, ne, np, layout) + offset
            ∂F∂x[ix] = get_partial(v, p_i)
        end
    end

    function store_partials!(∂F∂x::AbstractMatrix, v, i, ne, np, offset)
        for j in 1:size(∂F∂x, 2)
            for p_i in 1:np
                ix = alignment_linear_index(i, p_i, ne, np, layout) + offset
                ∂F∂x[ix, j] = get_partial(v[j], p_i)
            end
        end
    end

    function diff_entity!(∂F∂x, state, i, S, ne, np, offset)
        state_i = local_ad(state, i, S)
        v = F(eval_model, state_i, extra_arg...)
        store_partials!(∂F∂x, v, i, ne, np, offset)
    end

    offset = 0
    for e in get_primary_variable_ordered_entities(model)
        np = number_of_partials_per_entity(model, e)
        ne = count_active_entities(model.domain, e)
        if isnothing(sparsity)
            it_rng = 1:ne
        else
            it_rng = sparsity[e]
        end
        if length(it_rng) > 0
            ltag = get_entity_tag(tag, e)
            S = typeof(get_ad_entity_scalar(1.0, np, tag = ltag))
            for i in it_rng
                diff_entity!(∂F∂x, state, i, S, ne, np, offset)
            end
        end
        offset += ne*np
    end
end

function update_sensitivities!(∇G, i, G, adjoint_storage, state0, state, state_next, timesteps, all_forces)
    # Unpack simulators
    parameter_sim = adjoint_storage.parameter
    backward_sim = adjoint_storage.backward
    forward_sim = adjoint_storage.forward
    λ = adjoint_storage.lagrange
    λ_b = adjoint_storage.lagrange_buffer
    dparam = adjoint_storage.dparam
    # Timestep logic
    N = length(timesteps)
    forces = forces_for_timestep(forward_sim, all_forces, timesteps, i)
    dt = timesteps[i]
    # Assemble Jacobian w.r.t. current step
    t = sum(timesteps[1:i-1])
    @tic "jacobian (standard)" adjoint_reassemble!(forward_sim, state, state0, dt, forces, t)
    # Note the sign: There is an implicit negative sign in the linear solver when solving for the Newton increment. Therefore, the terms of the
    # right hand side are added with a positive sign instead of negative.
    rhs = adjoint_storage.rhs
    # Fill rhs with (∂J / ∂x)ᵀₙ (which will be treated with a negative sign when the result is written by the linear solver)
    S_p = get_objective_sparsity(adjoint_storage, :forward)
    @tic "objective primary gradient" state_gradient_outer!(rhs, G, forward_sim.model, forward_sim.storage.state, (dt, i, forces), sparsity = S_p)
    if isnothing(state_next)
        @assert i == N
        @. λ = 0
    else
        dt_next = timesteps[i+1]
        forces_next = forces_for_timestep(backward_sim, all_forces, timesteps, i+1)
        @tic "jacobian (with state0)" adjoint_reassemble!(backward_sim, state_next, state, dt_next, forces_next, t + dt_next)
        lsys_next = backward_sim.storage.LinearizedSystem
        op = linear_operator(lsys_next, skip_red = true)
        # In-place version of
        # rhs += op*λ
        # - (∂Fₙ₊₁ / ∂xₙ)ᵀ λₙ₊₁
        adjoint_transfer_canonical_order!(λ_b, λ, forward_sim.model, to_canonical = false)
        sens_add_mult!(rhs, op, λ_b)
    end
    # We have the right hand side, assemble the Jacobian and solve for the Lagrange multiplier
    lsolve = adjoint_storage.linear_solver
    dx = adjoint_storage.dx
    if adjoint_storage.multiple_rhs
        lsolve_arg = (dx = dx, r = rhs)
    else
        lsolve_arg = NamedTuple()
    end
    lsys = forward_sim.storage.LinearizedSystem
    if adjoint_storage.rhs_transfer_needed
        lsys.r_buffer .= rhs
    end
    @tic "linear solve" lstats = solve!(lsys, lsolve, forward_sim.model, forward_sim.storage; lsolve_arg...)
    adjoint_transfer_canonical_order!(λ, dx, forward_sim.model, to_canonical = true)
    # ∇ₚG = Σₙ (∂Fₙ / ∂p)ᵀ λₙ
    # Increment gradient
    @tic "jacobian (for parameters)" adjoint_reassemble!(parameter_sim, state, state0, dt, forces, t)
    lsys_param = parameter_sim.storage.LinearizedSystem
    op_p = linear_operator(lsys_param)
    sens_add_mult!(∇G, op_p, λ)

    @tic "objective parameter gradient" if !isnothing(dparam)
        if i == N
            @. dparam = 0
        end
        pbuf = adjoint_storage.param_buf
        state_gradient_outer!(pbuf, G, parameter_sim.model, parameter_sim.storage.state, (dt, i, forces))
        @. dparam += pbuf
    end
end

function sens_add_mult!(x::AbstractVector, op::LinearOperator, y::AbstractVector)
    mul!(x, op, y, 1.0, 1.0)
end

function sens_add_mult!(x::AbstractMatrix, op::LinearOperator, y::AbstractMatrix)
    for i in axes(x, 2)
        x_i = vec(view(x, :, i))
        y_i = vec(view(y, :, i))
        sens_add_mult!(x_i, op, y_i)
    end
end

function adjoint_reassemble!(sim, state, state0, dt, forces, time)
    s = sim.storage
    model = sim.model
    # Deal with state0 first
    reset_previous_state!(sim, state0)
    update_secondary_variables!(s, model, true)
    # Apply logic as if timestep is starting
    update_before_step!(s, model, dt, forces, time = time)
    # Then the current primary variables
    reset_variables!(s, model, state)
    update_state_dependents!(s, model, dt, forces)
    # Finally update the system
    update_linearized_system!(s, model)
end

function swap_primary_with_parameters!(pmodel, model, targets = parameter_targets(model))
    for (k, v) in pairs(model.primary_variables)
        if !(k in targets)
            set_parameters!(pmodel; k => v)
        end
    end
    # Original model holds the parameters, use those
    for (k, v) in pairs(model.parameters)
        if k in targets
            set_primary_variables!(pmodel; k => v)
        end
    end
    out = vcat(keys(model.secondary_variables)..., keys(model.primary_variables)...)
    for k in out
        push!(pmodel.output_variables, k)
    end
    unique!(pmodel.output_variables)
    return pmodel
end

parameter_targets(model::SimulationModel) = keys(get_parameters(model))

function adjoint_parameter_model(model, arg...; context = DefaultContext())
    # By default the adjoint model uses the default context since no linear solver
    # is needed.
    pmodel = adjoint_model_copy(model; context = context)
    # Swap parameters and primary variables
    swap_primary_with_parameters!(pmodel, model, arg...)
    return pmodel
end

function adjoint_model_copy(model::SimulationModel{O, S, F, C}; context = model.context) where {O, S, C, F}
    pvar = copy(model.primary_variables)
    svar = copy(model.secondary_variables)
    outputs = vcat(keys(pvar)..., keys(svar)...)
    prm = copy(model.parameters)
    extra = deepcopy(model.extra)
    eqs = model.equations
    # Transpose the system
    new_context = adjoint(context)
    N = typeof(new_context)
    return SimulationModel{O, S, F, N}(model.domain, model.system, new_context, model.formulation, model.plot_mesh, pvar, svar, prm, eqs, outputs, extra)
end

"""
    solve_numerical_sensitivities(model, states, reports, G, target;
                                                forces = setup_forces(model),
                                                state0 = setup_state(model),
                                                parameters = setup_parameters(model),
                                                epsilon = 1e-8)

Compute sensitivities of `model` parameter with name `target` for objective function `G`.

This method uses numerical perturbation and is primarily intended for testing of `solve_adjoint_sensitivities`.
"""
function solve_numerical_sensitivities(model, states, reports, G, target;
                                                forces = setup_forces(model),
                                                state0 = setup_state(model),
                                                parameters = setup_parameters(model),
                                                epsilon = 1e-8)
    timesteps = report_timesteps(reports)
    N = length(states)
    @assert length(reports) == N == length(timesteps)
    # Base objective
    base_obj = evaluate_objective(G, model, states, timesteps, forces)
    # Define perturbation
    param_var, param_num = get_parameter_pair(model, parameters, target)
    sz = size(param_num)
    grad_num = zeros(sz)
    grad_is_vector = grad_num isa AbstractVector
    if grad_is_vector
        grad_num = grad_num'
        sz = size(grad_num)
    end
    scale = variable_scale(param_var)
    if isnothing(scale)
        ϵ = epsilon
    else
        ϵ = scale*epsilon
    end
    n, m = sz
    for i in 1:n
        for j in 1:m
            param_i = deepcopy(parameters)
            perturb_parameter!(model, param_i, target, i, j, sz, ϵ)
            sim_i = Simulator(model, state0 = copy(state0), parameters = param_i)
            states_i, reports = simulate!(sim_i, timesteps, info_level = -1, forces = forces)
            v = evaluate_objective(G, model, states_i, timesteps, forces)
            grad_num[i, j] = (v - base_obj)/ϵ
        end
    end
    if grad_is_vector
        grad_num = grad_num'
    end
    return grad_num
end

function solve_numerical_sensitivities(model, states, reports, G; kwarg...)
    out = Dict()
    for k in keys(model.parameters)
        out[k] = solve_numerical_sensitivities(model, states, reports, G, k; kwarg...)
    end
end

function get_parameter_pair(model, parameters, target)
    return (model.parameters[target], parameters[target])
end

function perturb_parameter!(model, param_i, target, i, j, sz, ϵ)
    param_i[target][j + i*(sz[1]-1)] += ϵ
end

function evaluate_objective(G, model, states, timesteps, all_forces; large_value = 1e20)
    if length(states) < length(timesteps)
        # Failure: Put to a big value.
        @warn "Partial data passed, objective set to large value $large_value."
        obj = large_value
    else
        F = i -> G(model, states[i], timesteps[i], i, forces_for_timestep(nothing, all_forces, timesteps, i))
        obj = sum(F, eachindex(states))
    end
    return obj
end

function store_sensitivities(model, result, prm_map)
    out = Dict{Symbol, Any}()
    store_sensitivities!(out, model, result, prm_map)
    return out
end

function store_sensitivities!(out, model, result, prm_map; kwarg...)
    variables = model.primary_variables
    layout = matrix_layout(model.context)
    return store_sensitivities!(out, model, variables, result, prm_map, layout; kwarg...)
end

function store_sensitivities!(out, model, variables, result, prm_map, ::EquationMajorLayout)
    scalar_valued_objective = result isa AbstractVector
    if scalar_valued_objective
        N = 1
    else
        N = size(result, 2)
    end
    for (k, var) in pairs(variables)
        if !haskey(prm_map, k)
            continue
        end
        (; n, offset) = prm_map[k]
        m = degrees_of_freedom_per_entity(model, var)
        if n > 0
            rng = (offset+1):(offset+n)
            if scalar_valued_objective
                result_k = view(result, rng)
                v = extract_sensitivity_subset(result_k, var, n, m, offset)
            else
                # Grab each of the sensitivities and put them in an array for simplicity
                v = map(
                    i -> extract_sensitivity_subset(view(result, rng, i), var, n, m, offset + n*(i-1)),
                    1:N
                )
            end
            out[k] = v
        else
            out[k] = similar(result, 0)
        end
    end
    return out
end

function extract_sensitivity_subset(r, var, n, m, offset)
    if var isa ScalarVariable
        v = r
    else
        v = reshape(r, n ÷ m, m)'
    end
    v = collect(v)
    return v
end

gradient_vec_or_mat(n, ::Nothing) = zeros(n)
gradient_vec_or_mat(n, m) = zeros(n, m)

function variable_mapper(model::SimulationModel, type = :primary; targets = nothing, config = nothing, offset = 0)
    vars = get_variables_by_type(model, type)
    out = Dict{Symbol, Any}()
    for (t, var) in vars
        if !isnothing(targets) && !(t in targets)
            continue
        end
        var = vars[t]
        n = number_of_values(model, var)
        if !isnothing(config)
            lumping = config[t][:lumping]
            if !isnothing(lumping)
                N = maximum(lumping)
                @assert unique(lumping) == 1:N "Lumping for $t must include all values from 1 to $N"
                @assert length(lumping) == n "Lumping for $t must have one entry for each value if present"
                n = N
            end
        end
        out[t] = (n = n, offset = offset, scale = variable_scale(var))
        offset += n
    end
    return (out, offset)
end

function rescale_sensitivities!(dG, model, parameter_map)
    for (k, v) in parameter_map
        (; n, offset, scale) = v
        if !isnothing(scale)
            interval = (offset+1):(offset+n)
            if dG isa AbstractVector
                dG_k = view(dG, interval)
            else
                dG_k = view(dG, interval, :)
            end
            @. dG_k /= scale
        end
    end
end

function swap_variables(state, parameters, model; variables = true)
    p_keys = keys(get_parameters(model))
    v_keys = keys(get_variables(model))

    out = Dict{Symbol, Any}()
    # If variables: state should contain all fields, except those that are now in parameters
    if variables
        internal_swapper!(out, state, parameters, v_keys, p_keys)
    else
        internal_swapper!(out, parameters, state, p_keys, v_keys)
    end
    return out
end

function internal_swapper!(out, A, B, keep, skip)
    for (k, v) in A
        if !(k in skip)
            out[k] = v
        end
    end
    for (k, v) in B
        if k in keep
            out[k] = v
        end
    end
end

function adjoint_transfer_canonical_order!(λ, dx, model::MultiModel; to_canonical = true)
    offset = 0
    for (k, m) in pairs(model.models)
        ndof = number_of_degrees_of_freedom(m)
        ix = (offset+1):(offset+ndof)
        λ_i = view(λ, ix)
        dx_i = view(dx, ix)
        adjoint_transfer_canonical_order_inner!(λ_i, dx_i, m, matrix_layout(m.context), to_canonical)
        offset += ndof
    end
end

function adjoint_transfer_canonical_order!(λ, dx, model; to_canonical = true)
    adjoint_transfer_canonical_order_inner!(λ, dx, model, matrix_layout(model.context), to_canonical)
end

function adjoint_transfer_canonical_order_inner!(λ, dx, model, ::EquationMajorLayout, to_canonical)
    @. λ = dx
end

function adjoint_transfer_canonical_order_inner!(λ, dx, model, ::BlockMajorLayout, to_canonical)
    bz = 0
    for e in get_primary_variable_ordered_entities(model)
        if bz > 0
            error("Assumed that block major has a single entity group")
        end
        bz = degrees_of_freedom_per_entity(model, e)
    end
    n = length(dx) ÷ bz
    if to_canonical
        for b in 1:bz
            for i in 1:n
                λ[(b-1)*n + i] = dx[(i-1)*bz + b]
            end
        end
    else
        for b in 1:bz
            for i in 1:n
                λ[(i-1)*bz + b] = dx[(b-1)*n + i]
            end
        end
    end
end
