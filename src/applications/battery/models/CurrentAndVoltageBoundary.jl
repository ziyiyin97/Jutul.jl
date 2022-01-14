export CurrentAndVoltageSystem, CurrentAndVoltageDomain, CurrentForce, VoltageForce
export VoltageVAr, CurrentVar

struct CurrentAndVoltageSystem <: TervSystem end

struct CurrentAndVoltageDomain <: TervDomain end
active_entities(d::CurrentAndVoltageDomain, ::Any) = [1]

const CurrentAndVoltageModel  = SimulationModel{CurrentAndVoltageDomain,CurrentAndVoltageSystem}
function number_of_cells(::CurrentAndVoltageDomain) 1 end

# Driving force for the test equation
struct CurrentForce
    current
end

# Driving force for the test equation
struct VoltageForce
    current
end
#abstract type DiagonalEquation <: TervEquation end
# Equations
struct CurrentEquation <: DiagonalEquation
    equation
    function CurrentEquation(model, npartials::Integer; context = DefaultContext(), kwarg...)
        D = model.domain
        nc = number_of_cells(D)
        @assert nc == 1 # We use nc for clarity of the interface - but it should always be one!
        ne = 1 # Single, scalar equation
        e = CompactAutoDiffCache(ne, nc, npartials, context = context; kwarg...)
        new(e)
    end
end

# Equations
struct ControlEquation <: DiagonalEquation
    equation
    function ControlEquation(model, npartials::Integer; context = DefaultContext(), kwarg...)
        D = model.domain
        nc = number_of_cells(D)
        @assert nc == 1 # We use nc for clarity of the interface - but it should always be one!
        ne = 1 # Single, scalar equation
        e = CompactAutoDiffCache(ne, nc, npartials, context = context; kwarg...)
        new(e)
    end
end

function declare_sparsity(model, e::CurrentEquation, layout)
    return SparsePattern(1, 1, 1, 1, layout)
end

function declare_sparsity(model, e::ControlEquation, layout)
    return SparsePattern(1, 1, 1, 1, layout)
end

struct VoltVar <: ScalarVariable end
struct CurrentVar <: ScalarVariable end

function select_equations_system!(eqs, domain, system::CurrentAndVoltageSystem, formulation)
    eqs[:current_equation] = (CurrentEquation, 1)
    #eqs[:control_equation] = (ControlEquation, 1)
end



function update_equation!(eq::CurrentEquation, storage, model, dt)
 #   current = storage.state.Current
    phi = storage.state.Phi
    equation = get_entries(eq)
    @. equation = phi*1e-10
    #equation[:current_equation] = current
    #equation[:control_equation] = 0.0 
end

#function update_equation!(eq::ControlEquation, storage, model, dt)
#    current = storage.state.Current
#    phi = storage.state.Phi
#    equation = get_entries(eq)
#    @. equation = current
    #equation[:control_equation] = 0.0
#end

function build_forces(model::SimulationModel{G, S}; sources = nothing) where {G<:CurrentAndVoltageDomain, S<:CurrentAndVoltageSystem}
    return (sources = sources,)
end

function apply_forces_to_equation!(storage, model, eq::CurrentEquation, currentFun, time)
    #current = storage.state.Current
    equation = get_entries(eq)
    @. equation -= currentFun(time)
    #equation[:control_equations] = current-currentVal(time)
end

#function apply_forces_to_equation!(storage, model, eq::ControlEquation, voltageVal, time)
#    phi = storage.state.Phi
#    equation = get_entries(eq)
#    @. equation = phi - voltageVal(time)
#    #equation[:control_equiation] = phi - voltageVal(time)
#end

function update_cross_term!(ct::InjectiveCrossTerm, eq::CurrentEquation, target_storage, source_storage,# target_model, source_model, 
    target_model, 
    source_model, dt)
    X_T = target_storage.state.Phi#Voltage
    X_S = source_storage.state.Phi
    function f(X_S, X_T)
        (X_T - X_S)*1e-7
    end
    # Source term with AD context from source model - will end up as off-diagonal block
    @. ct.crossterm_source = f(X_S, value(X_T))
    # Source term with AD context from target model - will be inserted into equation
    @. ct.crossterm_target = f(value(X_S), X_T)
end

function update_cross_term!(ct::InjectiveCrossTerm, eq::Conservation{Charge}, target_storage, source_storage,# target_model, source_model, 
    target_model, 
    source_model, dt)
    X_T = target_storage.state.Phi
    X_S = source_storage.state.Phi#voltage
    function f(X_S, X_T)
        (X_T - X_S)*1e-7
    end
    # Source term with AD context from source model - will end up as off-diagonal block
    @. ct.crossterm_source = f(X_S, value(X_T))
    # Source term with AD context from target model - will be inserted into equation
    @. ct.crossterm_target = f(value(X_S), X_T)
end

function select_primary_variables!(S, domain, system::CurrentAndVoltageSystem, formulation)
    S[:Phi] = VoltVar()
    #S[:Current] = CurrentVar()
end

#function select_secondary_variables_system!(S, domain, system::CurrentAndVoltageSystem, formulation)
#end