export update_state_dependents!, check_convergence

function get_primary_variables(model::SimulationModel)
    return model.primary_variables
end

function get_secondary_variables(model::SimulationModel)
    return model.secondary_variables
end

function get_variables(model::SimulationModel)
    return merge(get_primary_variables(model), get_secondary_variables(model))
end

"""
Get only the entities where primary variables are present,
sorted by their order in the primary variables.
"""
function get_primary_variable_ordered_entities(model::SimulationModel)
    out = []
    current_entity = nothing
    for p in values(model.primary_variables)
        u = associated_entity(p)
        if u != current_entity
            # Note: We assume that primary variables for the same entity follows
            # each other. This is currently asserted in the model constructor.
            push!(out, u)
            current_entity = u
        end
    end
    return out
end

function number_of_partials_per_entity(model::SimulationModel, entity::TervUnit)
    n = 0
    for pvar in values(get_primary_variables(model))
        if associated_entity(pvar) == entity
            n += degrees_of_freedom_per_entity(model, pvar)
        end
    end
    return n
end

"""
Set up a state. You likely want to overload setup_state! instead of this one.
"""
function setup_state(model::TervModel, arg...)
    state = Dict{Symbol, Any}()
    setup_state!(state, model, arg...)
    return state
end

"""
Initialize primary variables and other state fields, given initial values as a Dict
"""
function setup_state!(state, model::TervModel, init_values::Dict = Dict())
    for (psym, pvar) in get_primary_variables(model)
        initialize_variable_value!(state, model, pvar, psym, init_values, need_value = true)
    end
    for (psym, svar) in get_secondary_variables(model)
        initialize_variable_value!(state, model, svar, psym, init_values, need_value = false)
    end
    initialize_extra_state_fields!(state, model)
end

"""
Add variables that need to be in state, but are never AD variables (e.g. phase status flag)
"""
function initialize_extra_state_fields!(state, model::TervModel)
    initialize_extra_state_fields_domain!(state, model, model.domain)
    initialize_extra_state_fields_system!(state, model, model.system)
    initialize_extra_state_fields_formulation!(state, model, model.formulation)
end

function initialize_extra_state_fields_domain!(state, model, domain)
    # Do nothing
end

function initialize_extra_state_fields_system!(state, model, system)
    # Do nothing
end

function initialize_extra_state_fields_formulation!(state, model, formulation)
    # Do nothing
end

"""
Allocate storage for the model. You should overload setup_storage! if you have a custom
definition.
"""
function setup_storage(model::TervModel; kwarg...)
    d = TervStorage()
    setup_storage!(d, model; kwarg...)
    return d
end

"""
Initialize the already allocated storage at the beginning of a simulation.
Use this to e.g. set up extra stuff in state0 needed for initializing the simulation loop.
"""
function initialize_storage!(storage, model::TervModel; initialize_state0 = true)
    if initialize_state0
        # Convert state and parameters to immutable for evaluation
        state0_eval = convert_to_immutable_storage(storage[:state0])
        param_eval = convert_to_immutable_storage(storage[:parameters])
        # Evaluate everything (with doubles) to make sure that possible 
        update_secondary_variables!(state0_eval, model, param_eval)
        # Create a new state0 with the desired/required outputs and
        # copy over those values before returning them back
        state0 = Dict()
        for key in model.output_variables
            v = state0_eval[key]
            if !isa(v, ConstantWrapper) && eltype(v)<:Real
                state0[key] = v
            end
        end
        storage[:state0] = state0
    end
    synchronize(model.context)
end

"""
Allocate storage for a given model. The storage consists of all dynamic quantities used in
the simulation. The default implementation allocates properties, equations and linearized system.
"""
function setup_storage!(storage, model::TervModel; setup_linearized_system = true,
                                                      state0 = setup_state(model),
                                                      parameters = setup_parameters(model),
                                                      tag = nothing,
                                                      kwarg...)
    if !isnothing(state0)
        storage[:parameters] = parameters
        storage[:state0] = state0
        storage[:state] = convert_state_ad(model, state0, tag)
        storage[:primary_variables] = reference_primary_variables(storage, model) 
    end
    setup_storage_model(storage, model)
    storage[:equations] = setup_equations(storage, model; tag = tag, kwarg...) 
    if setup_linearized_system
        storage[:LinearizedSystem] = setup_linearized_system!(storage, model)
        # We have the equations and the linearized system.
        # Give the equations a chance to figure out their place in the Jacobians.
        align_equations_to_linearized_system!(storage, model)
    end
end

function setup_storage_model(storage, model)
    setup_storage_domain!(storage, model, model.domain)
    setup_storage_system!(storage, model, model.system)
    setup_storage_formulation!(storage,  model, model.formulation)
end


function setup_storage_domain!(storage, model, domain)
    # Do nothing
end

function setup_storage_system!(storage, model, system)
    # Do nothing
end

function setup_storage_formulation!(storage,  model, formulation)
    # Do nothing
end

function reference_primary_variables(storage, model::TervModel; kwarg...)
    primaries = OrderedDict()
    state = storage[:state]
    for key in keys(get_primary_variables(model))
        primaries[key] = state[key]
    end
    return primaries
end

function setup_equations(storage, model::TervModel; kwarg...)
    # We use ordered dict since equation alignment with primary variables matter.
    eqs = OrderedDict()
    setup_equations!(eqs, storage, model; kwarg...)
    return eqs
end

function setup_equations!(eqs, storage, model::TervModel; tag = nothing, kwarg...)
    outstr = "Setting up $(length(model.equations)) groups of governing equations...\n"
    if !isnothing(tag)
        outstr = "$tag: "*outstr
    end
    counter = 1
    num_equations_total = 0
    for (sym, eq) in model.equations
        proto = eq[1]
        num = eq[2]
        if length(eq) > 2
            # We recieved extra kw-pairs to pass on
            extra = eq[3]
        else
            extra = []
        end
        e = proto(model, num; extra..., tag = tag, kwarg...)
        ne = number_of_entities(model, e)
        n = num*ne

        outstr *= "Group $counter/$(length(model.equations)) $(String(sym)) as $proto:\n\t → $num equations on each of $ne $(associated_entity(e)) for $n equations in total.\n"
        eqs[sym] = e
        counter += 1
        num_equations_total += n
    end
    outstr *= "$num_equations_total equations total distributed over $counter groups.\n"
    @debug outstr
end


function get_sparse_arguments(storage, model)
    layout = matrix_layout(model.context)
    return get_sparse_arguments(storage, model, layout)
end

function get_sparse_arguments(storage, model, layout::Union{EquationMajorLayout, UnitMajorLayout})
    ndof = number_of_degrees_of_freedom(model)
    eqs = storage[:equations]
    I = []
    J = []
    numrows = 0
    primary_entities = get_primary_variable_ordered_entities(model)
    for eq in values(eqs)
        numcols = 0
        for u in primary_entities
            S = declare_sparsity(model, eq, u, layout)
            if !isnothing(S)
                push!(I, S.I .+ numrows) # Row indices, offset by the size of preceeding equations
                push!(J, S.J .+ numcols) # Column indices, offset by the partials in entities we have passed
            end
            numcols += number_of_degrees_of_freedom(model, u)
        end
        @assert numcols == ndof
        # Number of equations correspond to number of rows
        numrows += number_of_equations(model, eq)
    end
    I = vcat(I...)
    J = vcat(J...)
    return SparsePattern(I, J, numrows, ndof, layout)
end

function get_sparse_arguments(storage, model, layout::BlockMajorLayout)
    eqs = storage[:equations]
    I = []
    J = []
    numrows = 0
    primary_entities = get_primary_variable_ordered_entities(model)
    block_size = degrees_of_freedom_per_entity(model, primary_entities[1])
    ndof = number_of_degrees_of_freedom(model) ÷ block_size
    for eq in values(eqs)
        numcols = 0
        eqs_per_entity = number_of_equations_per_entity(eq)
        for u in primary_entities
            dof_per_entity = degrees_of_freedom_per_entity(model, u)
            @assert dof_per_entity == eqs_per_entity == block_size "Block major layout only supported for square blocks."
            S = declare_sparsity(model, eq, u, layout)
            if !isnothing(S)
                push!(I, S.I .+ numrows) # Row indices, offset by the size of preceeding equations
                push!(J, S.J .+ numcols) # Column indices, offset by the partials in entities we have passed
            end
            numcols += count_active_entities(model.domain, u)
        end
        @assert numcols == ndof "Assumed square block, was $numcols x $ndof"
        numrows += number_of_entities(model, eq)
    end
    it = index_type(model.context)

    I = Vector{it}(vcat(I...))
    J = Vector{it}(vcat(J...))
    return SparsePattern(I, J, numrows, ndof, layout, block_size)
end

function get_sparse_arguments2(storage, model, layout::UnitMajorLayout)
    ndof = number_of_degrees_of_freedom(model)
    eqs = storage[:equations]
    I = []
    J = []
    numrows = 0
    numcols = 0
    primary_entities = get_primary_variable_ordered_entities(model)

    for (u_no, u) in enumerate(primary_entities)
        npartials = degrees_of_freedom_per_entity(model, u)
        nu = count_entities(model.domain, u)
        for (eq_no, eq) in enumerate(values(eqs))
            S = declare_sparsity(model, eq, u, layout)
            row_ix = (S.I-1)*u_no + 1
            col_ix = (S.J-1)*eq_no + 1
            if !isnothing(S)
                push!(I, row_ix .+ numrows)
                push!(J, col_ix .+ numcols)
            end
        end
        numrows += npartials*nu
        numcols += npartials*nu
    end
    I = vcat(I...)
    J = vcat(J...)
    return SparsePattern(I, J, numrows, ndof, layout)
end

function setup_linearized_system!(storage, model::TervModel)
    # Linearized system is going to have dimensions of
    # total number of equations x total number of primary variables
    if !haskey(storage, :equations)
        error("Unable to allocate linearized system - no equations found.")
    end
    # layout = matrix_layout(model.context)
    sparg = get_sparse_arguments(storage, model)
    lsys = setup_linearized_system(sparg, model)
    storage[:LinearizedSystem] = lsys
    # storage[:LinearizedSystem] = transfer(model.context, lsys)
    return lsys
end

function setup_linearized_system(sparse_arg, model)
    context = model.context
    LinearizedSystem(sparse_arg, context, matrix_layout(context))
end

function align_equations_to_linearized_system!(storage, model::TervModel; kwarg...)
    eqs = storage[:equations]
    jac = storage[:LinearizedSystem].jac
    align_equations_to_jacobian!(eqs, jac, model; kwarg...)
end

function align_equations_to_jacobian!(equations, jac, model; equation_offset = 0, variable_offset = 0)
    for key in keys(equations)
        eq = equations[key]
        align_to_jacobian!(eq, jac, model, equation_offset = equation_offset, variable_offset = variable_offset)
        equation_offset += number_of_equations(model, eq)
    end
    equation_offset
end

function allocate_array(context::TervContext, value, n...)
    tmp = transfer(context, value)
    return repeat(tmp, n...)
end

"""
Perform updates of everything that depends on the state.

This includes properties, governing equations and the linearized system
"""
function update_state_dependents!(storage, model::TervModel, dt, forces; time = NaN, update_secondary = true)
    if update_secondary
        @timeit "secondary variables" update_secondary_variables!(storage, model)
    end
    @timeit "equations" update_equations!(storage, model, dt)
    @timeit "forces" apply_forces!(storage, model, dt, forces; time = time)
    @timeit "boundary conditions" apply_boundary_conditions!(storage, model)
end

function apply_boundary_conditions!(storage, model::TervModel)
    parameters = storage.parameters
    apply_boundary_conditions!(storage, parameters, model)
end

apply_boundary_conditions!(storage, parameters, model) = nothing

function update_equations!(storage, model, dt = nothing)
    equations = storage.equations
    for key in keys(equations)
        @timeit "$key" update_equation!(equations[key], storage, model, dt)
    end
end

function update_linearized_system!(storage, model::TervModel; kwarg...)
    equations = storage.equations
    lsys = storage.LinearizedSystem
    update_linearized_system!(lsys, equations, model; kwarg...)
end

function update_linearized_system!(lsys, equations, model::TervModel; equation_offset = 0)
    r_buf = lsys.r_buffer
    for key in keys(equations)
        @timeit "$key" begin
            eq = equations[key]
            nz = lsys.jac_buffer
            N = number_of_equations(model, eq)
            n = number_of_equations_per_entity(eq)
            m = N ÷ n
            r = as_cell_major_matrix(r_buf, n, m, model, equation_offset)

            update_linearized_system_equation!(nz, r, model, eq)
            equation_offset += number_of_equations(model, eq)
        end
    end
end

function check_convergence(storage, model; kwarg...)
    lsys = storage.LinearizedSystem
    eqs = storage.equations
    check_convergence(lsys, eqs, storage, model; kwarg...)
end

function check_convergence(lsys, eqs, storage, model; iteration = nothing, extra_out = false, tol = nothing, offset = 0, kwarg...)
    converged = true
    e = 0
    eoffset = 0
    r_buf = lsys.r_buffer
    prm = storage.parameters.tolerances
    if isnothing(tol)
        tol = prm.default
    end
    output = []
    for key in keys(eqs)
        eq = eqs[key]
        N = number_of_equations(model, eq)
        n = number_of_equations_per_entity(eq)
        m = N ÷ n
        r_v = as_cell_major_matrix(r_buf, n, m, model, offset)

        @timeit "$key" all_crits = convergence_criterion(model, storage, eq, r_v; kwarg...)
        e_keys = keys(all_crits)
        tols = Dict()
        for e_k in e_keys
            if haskey(prm, key)
                v = prm[key]
                if isa(v, AbstractFloat)
                    t_e = v
                else
                    t_e = v[e_k]
                end
            else
                t_e = tol
            end
            errors = all_crits[e_k].errors
            ok = errors .< t_e
            converged = converged && all(ok)
            e = maximum([e, maximum(errors)/t_e])
            tols[e_k] = t_e
        end
        offset += N
        eoffset += n
        if extra_out
            push!(output, (name = key, criterions = all_crits, tolerances = tols))
        end
    end
    if extra_out
        return (converged, e, output)
    else
        return converged
    end
end

"""
Apply a set of forces to all equations. Equations that don't support a given force
will just ignore them, thanks to the power of multiple dispatch.
"""
function apply_forces!(storage, model, dt, forces; time = NaN)
    equations = storage.equations
    for key in keys(equations)
        eq = equations[key]
        for fkey in keys(forces)
            force = forces[fkey]
            apply_forces_to_equation!(storage, model, eq, force, time)
        end
    end
end

function apply_forces!(storage, model, dt, ::Nothing; time = NaN)

end

function setup_parameters(model)
    d = Dict{Symbol, Any}()
    d[:tolerances] = Dict{Symbol, Any}()
    d[:tolerances][:default] = 1e-3
    setup_parameters_domain!(d, model, model.domain)
    setup_parameters_system!(d, model, model.system)
    setup_parameters_context!(d, model, model.context)
    setup_parameters_formulation!(d, model, model.formulation)
    return d
end
setup_parameters_domain!(d, model, ::Any) = nothing
setup_parameters_system!(d, model, ::Any) = nothing
setup_parameters_context!(d, model, ::Any) = nothing
setup_parameters_formulation!(d, model, ::Any) = nothing

function build_forces(model::TervModel)
    return NamedTuple()
end

function solve_and_update!(storage, model::TervModel, dt = nothing; linear_solver = nothing, recorder = nothing, kwarg...)
    lsys = storage.LinearizedSystem
    t_solve = @elapsed @timeit "linear solve" solve!(lsys, linear_solver, model, storage, dt, recorder)
    t_update = @elapsed @timeit "primary variables" update_primary_variables!(storage, model; kwarg...)
    return (t_solve, t_update)
end

function update_primary_variables!(storage, model::TervModel; kwarg...)
    dx = storage.LinearizedSystem.dx_buffer
    update_primary_variables!(storage.primary_variables, dx, model; kwarg...)
end

function update_primary_variables!(primary_storage, dx, model::TervModel; check = false)
    layout = matrix_layout(model.context)
    cell_major = is_cell_major(layout)
    offset = 0
    primary = get_primary_variables(model)
    ok = true
    if cell_major
        offset = 0 # Offset into global r array
        for u in get_primary_variable_ordered_entities(model)
            np = number_of_partials_per_entity(model, u)
            nu = count_entities(model.domain, u)
            t = isa(layout, BlockMajorLayout)
            if t
                Dx = get_matrix_view(dx, np, nu, true, offset)
            else
                Dx = get_matrix_view(dx, np, nu, false, offset)'
            end
            local_offset = 0
            for (pkey, p) in primary
                # This is a bit inefficient
                if u != associated_entity(p)
                    continue
                end
                ni = degrees_of_freedom_per_entity(model, p)
                dxi = view(Dx, :, (local_offset+1):(local_offset+ni))
                if check
                    ok_i = check_increment(dxi, p, pkey)
                    ok = ok && ok_i
                end
                @timeit "$pkey" update_primary_variable!(primary_storage, p, pkey, model, dxi)
                local_offset += ni
            end
            offset += nu*np
        end
    else
        for (pkey, p) in primary
            n = number_of_degrees_of_freedom(model, p)
            rng = (1:n) .+ offset
            dxi = view(dx, rng)
            if check
                ok_i = check_increment(dxi, p, pkey)
                ok = ok && ok_i
            end
            @timeit "$pkey" update_primary_variable!(primary_storage, p, pkey, model, dxi)
            offset += n
        end
    end
    if !ok
        error("Primary variables recieved invalid updates.")
    end
end

"""

"""
function update_before_step!(storage, model, dt, forces)
    update_before_step_domain!(storage, model, model.domain, dt, forces)
    update_before_step_system!(storage, model, model.system, dt, forces)
    update_before_step_formulation!(storage, model, model.formulation, dt, forces)
end

function update_before_step_domain!(state, model, domain, dt, forces)
    # Do nothing
end

function update_before_step_system!(state, model, system, dt, forces)
    # Do nothing
end

function update_before_step_formulation!(state, model, formulation, dt, forces)
    # Do nothing
end

function update_after_step!(storage, model, dt, forces)
    state = storage.state
    state0 = storage.state0
    for key in keys(state0)
        @. state0[key] = value(state[key])
    end
end

function get_output_state(storage, model)
    # As this point (after a converged step) state0 should be state without AD.
    s0 = storage.state0
    D = Dict{Symbol, Any}()
    for k in keys(s0)
        D[k] = copy(s0[k])
    end
    return D
end

function reset_to_previous_state!(storage, model)
    primary = storage.primary_variables
    state0 = storage.state0
    for f in keys(primary)
        if haskey(state0, f)
            update_values!(primary[f], state0[f])
        end
    end
end
