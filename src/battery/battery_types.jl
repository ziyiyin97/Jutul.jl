using Terv

# TODO: These are not all needed
export ElectroChemicalComponent, CurrentCollector, Electectrolyte, TestElyte
export vonNeumannBC, DirichletBC, BoundaryCondition, MinimalECTPFAGrid
export ChargeFlow, MixedFlow, Conservation, BoundaryPotential, BoundaryCurrent
export Phi, C, T, ChargeAcc, MassAcc, EnergyAcc, KGrad
export BOUNDARY_CURRENT, corr_type

###########
# Classes #
###########

abstract type ElectroChemicalComponent <: TervSystem end
# Alias for a genereal Electro Chemical Model
const ECModel = SimulationModel{<:Any, <:ElectroChemicalComponent, <:Any, <:Any}

abstract type ElectroChemicalGrid <: TervGrid end

# Potentials
# ? introduce abstract type Potential?
struct Phi <: ScalarVariable end
struct C <: ScalarVariable end
struct T <: ScalarVariable end

struct Conductivity <: ScalarVariable end
struct Diffusivity <: ScalarVariable end
struct ThermalConductivity <: ScalarVariable end

struct Conservation{T} <: TervEquation
    accumulation::TervAutoDiffCache
    accumulation_symbol::Symbol
    half_face_flux_cells::TervAutoDiffCache
    density::TervAutoDiffCache
    flow_discretization::FlowDiscretization
end

# Accumulation variables
abstract type AccumulationVariable <: ScalarVariable end
struct ChargeAcc <: AccumulationVariable end
struct MassAcc <: AccumulationVariable end
struct EnergyAcc <: AccumulationVariable end

# Currents corresponding to a accumulation type
const BOUNDARY_CURRENT = Dict(
    ChargeAcc() => :BCCharge,
    MassAcc()   => :BCMass,
    EnergyAcc() => :BCEnergy,
)
# TODO: Can this not be automated????
function corr_type(::Conservation{ChargeAcc}) ChargeAcc() end
function corr_type(::Conservation{MassAcc}) MassAcc() end
function corr_type(::Conservation{EnergyAcc}) EnergyAcc() end

# Represents k∇T, where k is a tensor, T a potential
abstract type KGrad{T} <: ScalarVariable end
struct TPkGrad{T} <: KGrad{T} end

abstract type ECFlow <: FlowType end
struct ChargeFlow <: ECFlow end
struct MixedFlow <: ECFlow end

struct BoundaryPotential{T} <: ScalarVariable end
struct BoundaryCurrent{T} <: ScalarVariable 
    cells
end

abstract type Current <: ScalarVariable end
struct TotalCurrent <: Current end
struct ChargeCarrierFlux <: Current end
struct EnergyFlux <: Current end

function number_of_units(model, pv::Current)
    return 2*count_units(model.domain, Faces())
end

abstract type NonDiagCellVariables <: TervVariables end

# Abstract type of a vector that is defined on a cell, from face flux
abstract type CellVector <: NonDiagCellVariables end
struct JCell <: CellVector end

abstract type ScalarNonDiagVaraible <: NonDiagCellVariables end
struct JSq <: ScalarNonDiagVaraible end

struct JSqDiag <: ScalarVariable end

struct MinimalECTPFAGrid{R<:AbstractFloat, I<:Integer} <: ElectroChemicalGrid
    """
    Simple grid for a electro chemical component
    """
    volumes::AbstractVector{R}
    neighborship::AbstractArray{I}
    boundary_cells::AbstractArray{I}
    boundary_T_hf::AbstractArray{R}
    P::AbstractArray{R} # Tensor to map from cells to faces
    S::AbstractArray{R} # Tensor map cell vector to cell scalar
    cellcellvectbl
    cellcelltbl
end

################
# Constructors #
################


function MinimalECTPFAGrid(pv, N, bc=[], T_hf=[], P=[], S=[])
    nc = length(pv)
    pv::AbstractVector
    @assert size(N, 1) == 2
    if length(N) > 0
        @assert minimum(N) > 0
        @assert maximum(N) <= nc
    end
    @assert all(pv .> 0)
    @assert size(bc) == size(T_hf)
    
    cellcellvectbl = get_cellcellvec_map(N)
    cellcelltbl  = get_cellcell_map(N)

    MinimalECTPFAGrid{eltype(pv), eltype(N)}(
        pv, N, bc, T_hf, P, S, cellcellvectbl, cellcelltbl
        )
end

function Conservation(
    acc_type, model, number_of_equations;
    flow_discretization = model.domain.discretizations[1], kwarg...
    )
    """
    A conservation law corresponding to the underlying potential pvar
    """

    accumulation_symbol = acc_symbol(acc_type)

    D = model.domain
    cell_unit = Cells()
    face_unit = Faces()
    nc = count_units(D, cell_unit)
    nf = count_units(D, face_unit)
    nhf = 2 * nf
    nn = size(get_neighborship(D.grid), 2)
    n_tot = nc + 2 * nn

    alloc = (n, unit, n_units_pos) -> CompactAutoDiffCache(
        number_of_equations, n, model, unit = unit, n_units_pos = n_units_pos,
        context = model.context; kwarg...
    )

    acc = alloc(nc, cell_unit, nc)
    hf_cells = alloc(nhf, cell_unit, nhf)
    density = alloc(n_tot, cell_unit, n_tot)

    Conservation{typeof(acc_type)}(
        acc, accumulation_symbol, hf_cells, density, flow_discretization
    )
end


function acc_symbol(::ChargeAcc)
    return :ChargeAcc
end

function acc_symbol(::MassAcc)
    return :MassAcc
end

function acc_symbol(::EnergyAcc)
    return :EnergyAcc
end
