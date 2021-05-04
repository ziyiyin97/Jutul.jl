export MultiPhaseSystem, ImmiscibleSystem, SinglePhaseSystem
export LiquidPhase, VaporPhase
export number_of_phases, get_short_name, get_name, subscript
export update_linearized_system!
export SourceTerm
export setup_state, setup_state!

export allocate_storage, update_equations!

using CUDA
# Abstract multiphase system
abstract type MultiPhaseSystem <: TervSystem end


function get_phases(sys::MultiPhaseSystem)
    return sys.phases
end


function number_of_phases(sys::MultiPhaseSystem)
    return length(get_phases(sys))
end

struct SourceTerm{R<:Real,I<:Integer}
    cell::I
    values::AbstractVector{R}
end


## Systems
# Immiscible multiphase system
struct ImmiscibleSystem <: MultiPhaseSystem
    phases::AbstractVector
end

# function ImmiscibleSystem(phases)
#    @assert length(phases) > 1 "System should have at least two phases. For single-phase, use SinglePhaseSystem instead."
# end

# Single-phase
struct SinglePhaseSystem <: MultiPhaseSystem
    phase
end

function get_phases(sys::SinglePhaseSystem)
    return [sys.phase]
end

function number_of_phases(::SinglePhaseSystem)
    return 1
end

## Phases
# Abstract phase
abstract type AbstractPhase end

function get_short_name(phase::AbstractPhase)
    return get_name(phase)[1:1]
end

function subscript(prefix::String, phase::AbstractPhase)
    return string(prefix, "_", get_short_name(phase))
end
# Liquid phase
struct LiquidPhase <: AbstractPhase end

function get_name(::LiquidPhase)
    return "Liquid"
end

# Vapor phases
struct VaporPhase <: AbstractPhase end

function get_name(::VaporPhase)
    return "Vapor"
end

## Main implementation
function setup_state(model, arg...)
    state = Dict{String, Any}()
    setup_state!(state, model, model.grid, model.system, arg...)
    return state
end

function setup_state!(state, model, G, sys::MultiPhaseSystem, init_values)
    pvars = get_primary_variables(model)
    for pvar in get_primary_variables(model)
        initialize_primary_variable_value(state, model, pvar, init_values)
    end

    nc = number_of_cells(model.grid)
    phases = get_phases(sys)
    for phase in phases
        state[subscript("PhaseMass", phase)] = transfer(model.context, zeros(nc))
    end
    if !isa(sys, SinglePhaseSystem)
        s = state["Saturations"]
        for (i, phase) in enumerate(phases)
            state[subscript("Saturation", phase)] = view(s, i, :)
        end
    end
end

function update_state!(model, storage)
    lsys = storage["LinearizedSystem"]
    state = storage["state"]

    offset = 0
    primary = get_primary_variables(model)
    for p in primary
        n = number_of_degrees_of_freedom(model, p)
        rng = (1:n) .+ offset
        update_state!(state, p, model, view(lsys.dx, rng))
        offset += n
    end
end


function convert_state_ad(model, state)
    context = model.context
    stateAD = deepcopy(state)
    vars = String.(keys(state))

    primary = get_primary_variables(model)
    # Loop over primary variables and set them to AD, with ones at the correct diagonal
    counts = map((x) -> degrees_of_freedom_per_unit(x), primary)
    n_partials = sum(counts)
    @debug "Found $n_partials primary variables."
    offset = 0
    for (i, pvar) in enumerate(primary)
        stateAD = initialize_primary_variable_ad(stateAD, model, pvar, offset, n_partials)
        offset += counts[i]
    end
    primary_names = get_primary_variable_names(model)
    secondary = setdiff(vars, primary_names)
    # Loop over secondary variables and initialize as AD with zero partials
    for s in secondary
        stateAD[s] = allocate_array_ad(stateAD[s], context = context, npartials = n_partials)
    end
    return stateAD
end

# Primary variable logic

# Pressure as primary variable
struct Pressure <: ScalarPrimaryVariable
    name
end

function Pressure()
    Pressure("Pressure")
end
# Saturations as primary variable
struct Saturations <: GroupedPrimaryVariables
    name
    phases
end

function Saturations(phases::AbstractArray)
    Saturations("Saturations", phases)
end

function initialize_primary_variable_ad(state, model, pvar::Saturations, offset, npartials)
    name = get_name(pvar)
    nph = length(pvar.phases)
    # nph - 1 primary variables, with the last saturation being initially zero AD
    dp = vcat((1:nph-1) .+ offset, 0)
    v = state[name]
    v = allocate_array_ad(v, diag_pos = dp, context = model.context, npartials = npartials)
    for i in 1:size(v, 2)
        v[end, i] = 1 - sum(v[1:end-1, i])
    end
    state[name] = v
    return state
end

function initialize_primary_variable_value(state, model, pvar::Saturations, val)
    name = get_name(pvar)
    if isa(val, Dict)
        val = val[name]
    end
    val::AbstractVecOrMat # Should be a vector or a matrix
    V = deepcopy(val)
    @assert size(V, 1) == number_of_phases(model.system)
    n = number_of_units(model, pvar)
    if isa(V, AbstractVector)
        V = repeat(V, 1, n)
    end
    @assert size(V, 2) == n
    state[name] = transfer(model.context, V)
    return state
end

function degrees_of_freedom_per_unit(v::Saturations)
    return length(v.phases) - 1
end

function get_names(v::Saturations)
    return map((x) -> subscript(v.name, x), v.phases[1:end-1])
end

function get_name(v::Saturations)
    return v.name
end

function select_primary_variables(system::SinglePhaseSystem, formulation, discretization)
    return [Pressure()]
end

function select_primary_variables(system::ImmiscibleSystem, formulation, discretization)
    return [Pressure(), Saturations(system.phases)]
end

function allocate_storage!(d, model::SimulationModel{T, S}) where {T<:Any, S<:MultiPhaseSystem}
    G = model.grid
    sys = model.system
    context = model.context

    nph = number_of_phases(sys)
    phases = get_phases(sys)
    npartials = nph
    nc = number_of_cells(G)
    nhf = number_of_half_faces(G)

    jac = get_sparsity_pattern(G, nph, nph)

    n_dof = nc*nph
    dx = zeros(n_dof)
    r = zeros(n_dof)
    lsys = LinearizedSystem(jac, r, dx)
    alloc = (n) -> allocate_array_ad(nph, n, context = context, npartials = npartials)
    alloc_value = (n) -> allocate_array_ad(nph, n, context = context, npartials = 0)
    law = ConservationLaw(G, lsys, npartials, context = context, equations_per_unit = nph)
    d["MassConservation"] = law
    # Mobility of phase
    d["Mobility"] = alloc(nc)
    # Mass density of phase
    d["Density"] = alloc(nc)
    # Mobility * Density. We compute store this separately since density
    # and mobility are both used in other spots
    d["MassMobility"] = alloc(nc)
    d["TotalMass"] = alloc(nc)
    d["TotalMass0"] = alloc_value(nc)

    # Transfer linearized system afterwards since the above manipulations are
    # easier to do on CPU
    d["LinearizedSystem"] = transfer(context, lsys)
end

function initialize_storage!(d, model::SimulationModel{T, S}) where {T<:Any, S<:MultiPhaseSystem}
    update_properties!(model, d)
    m = d["TotalMass"]
    m0 = d["TotalMass0"]
    @. m0 = value(m)
end

function update_equations!(model::SimulationModel{G, S}, storage; 
    dt = nothing, sources = nothing) where {G<:MinimalTPFAGrid, S<:MultiPhaseSystem}
    update_properties!(model, storage)
    acc = update_accumulation!(model, storage, dt)
    phases = get_phases(model.system)
    for phNo in eachindex(phases)
        if !isnothing(sources)
            # @debug "Inserting source terms."
            insert_sources(acc, sources, phNo)
        end
    end
    update_half_face_flux!(model, storage)
end

function update_properties!(model, storage)
    update_density!(model, storage)
    update_mobility!(model, storage)
    update_mass_mobility!(model, storage)
    update_total_mass!(model, storage)
end

# Updates of various cell properties follows
function update_mobility!(model::SimulationModel{G, S}, storage) where {G<:Any, S<:SinglePhaseSystem}
    phase = model.system.phase
    mob = storage["Mobility"]
    mu = storage["parameters"][subscript("Viscosity", phase)]
    fapply!(mob, () -> 1/mu)
end

function update_mobility!(model::SimulationModel{G, S}, storage) where {G<:Any, S<:ImmiscibleSystem}
    p = storage["parameters"]
    mob = storage["Mobility"]
    n = p["CoreyExponents"]
    mu = p["Viscosity"]

    s = storage["state"]["Saturations"]
    @. mob = s^n/mu
end

function update_density!(model, storage)
    param = storage["parameters"]
    state = storage["state"]
    p = state["Pressure"]
    rho = storage["Density"]
    phases = get_phases(model.system)
    for (i, phase) in enumerate(phases)
        d = subscript("Density", phase)
        rho_i = view(rho, i, :)
        # rho = storage[d]
        r = param[d]
        if isa(r, NamedTuple)
            f_rho = (p) -> r.rhoS*exp((p - r.pRef)*r.c)
        else
            # Function handle
            f_rho = r
        end
        fapply!(rho_i, f_rho, p)
    end
    return rho
end

function update_total_mass!(model::SimulationModel{G, S}, storage) where {G<:Any, S<:SinglePhaseSystem}
    pv = model.grid.pv
    rho = storage["Density"]
    totMass = storage["TotalMass"]
    fapply!(totMass, *, rho, pv)
end

function update_total_mass!(model::SimulationModel{G, S}, storage) where {G<:Any, S<:ImmiscibleSystem}
    pv = model.grid.pv
    rho = storage["Density"]
    totMass = storage["TotalMass"]
    s = storage["state"]["Saturations"]
    fapply!(totMass, *, rho, pv, s)
end


function update_total_mass!(model::SimulationModel{G, S}, storage, phase::AbstractPhase) where {G<:Any, S<:ImmiscibleSystem}
    pv = model.grid.pv
    rho = storage[subscript("Density", phase)]
    s = storage["state"][subscript("Saturation", phase)]
    totMass = storage[subscript("TotalMass", phase)]
    fapply!(totMass, *, rho, pv, s)
end

function update_mass_mobility!(model, storage)
    mobrho = storage["MassMobility"]
    mob = storage["Mobility"]
    rho = storage["Density"]
    # Assign the values
    fapply!(mobrho, *, mob, rho)
end

# Update of discretization terms
function update_accumulation!(model, storage, dt)
    law = storage["MassConservation"]
    acc = law.accumulation
    mass = storage["TotalMass"]
    mass0 = storage["TotalMass0"]
    fapply!(acc, (m, m0) -> (m - m0)/dt, mass, mass0)
    return acc
end

function update_half_face_flux!(model, storage)
    law = storage["MassConservation"]
    p = storage["state"]["Pressure"]
    mmob = storage["MassMobility"]

    half_face_flux!(model, law.half_face_flux, mmob, p)
end

# Source terms, etc
@inline function insert_sources(acc, sources, phNo)
    for src in sources
        @inbounds acc[src.cell] -= src.values[phNo]
    end
end

@inline function insert_sources(acc::CuArray, sources, phNo)
    s = cu(map(x -> x.values[phNo], sources))
    i = cu(map(x -> x.cell, sources))
    @. acc[i] -= s
end

# Updating of linearized system
function update_linearized_system!(model::TervModel, storage)
    sys = model.system;
    sys::MultiPhaseSystem

    lsys = storage["LinearizedSystem"]
    update_linearized_system!(model, lsys, storage["MassConservation"])
end

function update_linearized_system!(model, lsys::LinearizedSystem, law::ConservationLaw)
    G = model.grid
    context = model.context
    ker_compat = kernel_compatibility(context)
    apos = law.accumulation_jac_pos
    jac = lsys.jac
    r = lsys.r
    # Fill in diagonal
    # @info "Accumulation fillin"
    fill_accumulation!(jac, r, law.accumulation, apos, context, ker_compat)
    # Fill in off-diagonal
    fpos = law.half_face_flux_jac_pos
    # @info "Half face flux fillin"
    fill_half_face_fluxes(jac, r, G.conn_pos, law.half_face_flux, apos, fpos, context, ker_compat)
end

# Accumulation: Base implementation
"Fill acculation term onto diagonal with pre-determined access pattern into jac"
function fill_accumulation!(jac, r, acc, apos, context, ::KernelDisallowed)
    nzval = get_nzval(jac)
    nder = size(apos, 1)
    @inbounds Threads.@threads for col = 1:size(apos, 2)
        r[col] = acc[col].value
        fill_accumulation_jac!(nzval, acc, apos, col, nder)
    end
end

@inline function fill_accumulation_jac!(nzval, acc, apos, col, nder)
    @inbounds for derNo = 1:nder
        nzval[apos[derNo, col]] = acc[col].partials[derNo]
    end
end
# Kernel / CUDA version follows
function fill_accumulation!(jac, r, acc, apos, context, ::KernelAllowed)
    @. r = value(acc)
    jz = get_nzval(jac)
    @inbounds for i = 1:size(apos, 1)
        jz[apos[i, :]] = map((x) -> x.partials[i], acc)
    end
    # kernel = fill_accumulation_jac_kernel(context.device, context.block_size)
    # event = kernel(jz, acc, apos, ndrange = size(apos))
    # wait(event)
end

"Kernel for filling Jacobian from accumulation term"
@kernel function fill_accumulation_jac_kernel(nzval, @Const(acc), @Const(apos))
    derNo, col = @index(Global, NTuple)
    @inbounds nzval[apos[derNo, col]] = acc[col].partials[derNo]
    # i = @index(Global, Linear)
    # @inbounds nzval[apos[1, i]] = acc[i].partials[1]
end

# Fluxes: Base implementation
"Fill fluxes onto diagonal with pre-determined access pattern into jac. Essentially performs Div ( flux )"
function fill_half_face_fluxes(jac, r, conn_pos, half_face_flux, apos, fpos, context, ::KernelDisallowed)
    Jz = get_nzval(jac)
    Threads.@threads for col = 1:length(apos)
        @inbounds for i = conn_pos[col]:(conn_pos[col+1]-1)
            # Update diagonal value
            r[col] += half_face_flux[i].value
            # Fill inn Jacobian values
            fill_half_face_fluxes_jac!(Jz, half_face_flux, apos, fpos, col, i)
        end
    end
end

@inline function fill_half_face_fluxes_jac!(nzval, half_face_flux, apos, fpos, col, i)
    @inbounds for derNo = 1:size(apos, 1)
        index = fpos[derNo, i]
        diag_index = apos[derNo, col]
        df_di = half_face_flux[i].partials[derNo]
        nzval[index] = -df_di
        nzval[diag_index] += df_di
    end
end

# Kernel / CUDA version follows
function fill_half_face_fluxes(jac, r, conn_pos, half_face_flux, apos, fpos, context, ::KernelAllowed)
    d = size(apos)
    Jz = get_nzval(jac)

    # println(":: HF, value")
    begin
        ncol = d[2]
        kernel = fill_half_face_fluxes_val_kernel(context.device, context.block_size, ncol)
        event_val = kernel(r, half_face_flux, conn_pos, ndrange = ncol)
        # wait(event_val)
    end

    # println(":: HF, offdiag")
    if true
        begin
            @inbounds for i = 1:size(fpos, 1)
                Jz[fpos[i, :]] = map((x) -> -x.partials[i], half_face_flux)
            end
        end
        # println(":: HF, diag")
        begin
            kernel = fill_half_face_fluxes_jac_diag_kernel(context.device, context.block_size, d)
            event_jac = kernel(Jz, half_face_flux, apos, conn_pos, ndrange = d)
            # wait(event_jac)
        end
        wait(event_val)
        wait(event_jac)
    else
        println(":: HF, full kernel")
        @time begin
            kernel = fill_half_face_fluxes_jac_kernel(context.device, context.block_size, d)
            event_jac = kernel(Jz, half_face_flux, apos, fpos, conn_pos, ndrange = d)
            wait(event_jac)
        end
    end
end

@kernel function fill_half_face_fluxes_val_kernel(r, @Const(half_face_flux), @Const(conn_pos))
    col = @index(Global, Linear)
    v = 0
    @inbounds for i = conn_pos[col]:(conn_pos[col+1]-1)
        # Update diagonal value
        v += half_face_flux[i].value
    end
    @inbounds r[col] += v
end


@kernel function fill_half_face_fluxes_jac_diag_kernel(nzval, @Const(half_face_flux), @Const(apos), @Const(conn_pos))
    derNo, col = @index(Global, NTuple)
    v = 0
    @inbounds for i = conn_pos[col]:(conn_pos[col+1]-1)
        v += half_face_flux[i].partials[derNo]
    end
    @inbounds nzval[apos[derNo, col]] += v
    nothing
end

@kernel function fill_half_face_fluxes_jac_kernel(nzval, @Const(half_face_flux), @Const(apos), @Const(fpos), @Const(conn_pos))
    derNo, col = @index(Global, NTuple)
    v = 0
    @inbounds for i = conn_pos[col]:(conn_pos[col+1]-1)
        index = fpos[derNo, i]
        df_di = half_face_flux[i].partials[derNo]
        nzval[index] = -df_di
        v += df_di
    end
    @inbounds diag_index = apos[derNo, col]
    @inbounds nzval[diag_index] += v # Should this be atomic?
    nothing
end

@inline function get_nzval(jac)
    return jac.nzval
end

@inline function get_nzval(jac::AbstractCuSparseMatrix)
    # Why does CUDA and Base differ on capitalization?
    return jac.nzVal
end
