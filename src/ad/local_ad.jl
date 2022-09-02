import Base: getindex, @propagate_inbounds, parent, size, axes
using ForwardDiff

struct LocalPerspectiveAD{T, N, A<:AbstractArray{T,N}, I} <: AbstractArray{T,N}
    index::I
    data::A
end

function LocalPerspectiveAD(a::A, index::I_t) where {A<:AbstractArray, I_t<:Integer}
    LocalPerspectiveAD{eltype(a), ndims(A), A, I_t}(index, a)
end

struct LocalStateAD{T, I, E} # Data type, index, entity tag
    index::I
    data::T
end

struct ValueStateAD{T} # Data type
    data::T
end

as_value(x::Union{NamedTuple,AbstractDict,JutulStorage}) = ValueStateAD(x)

export local_ad
@inline local_ad(v::AbstractArray, i::Integer) = LocalPerspectiveAD(v, i)
@inline local_ad(v::ConstantWrapper, i::Integer) = v
@inline local_ad(v, ::Nothing) = as_value(v)
@inline local_ad(v, i) = v


@inline function new_entity_index(state::LocalStateAD{T, I, E}, index::I) where {T, I, E}
    return LocalStateAD{T, I, E}(index, getfield(state, :data))
end

@inline function new_entity_index(x, index)
    return x
end


@inline local_entity(a::LocalPerspectiveAD) = a.index

@inline function value_or_ad(A::LocalPerspectiveAD{T}, v::T, entity) where T
    if entity === local_entity(A)
        return v
    else
        return T(value(v))
    end
end

@inline @propagate_inbounds function Base.getindex(A::LocalPerspectiveAD{T}, i::Int) where T
    d = A.data[i]
    return value_or_ad(A, d, i)
end

@inline @propagate_inbounds function Base.getindex(A::LocalPerspectiveAD{T}, i::Int, j::Int) where T
    d = A.data[i, j]
    return value_or_ad(A, d, j)
end

@inline Base.parent(A::LocalPerspectiveAD) = A.data
@inline Base.size(A::LocalPerspectiveAD) = size(A.data)
@inline Base.axes(A::LocalPerspectiveAD) = axes(A.data)
@inline parenttype(::Type{LocalPerspectiveAD{T,N,A,I}}) where {T,N,A,I} = A
@inline Base.haskey(state::LocalStateAD, f::Symbol) = haskey(getfield(state, :data), f)


# Match in type - pass index on
@inline next_level_local_ad(x::AbstractArray{T}, ::Type{T}, index) where T = local_ad(x, index)
# Mismatch in AD type - take value
@inline next_level_local_ad(x, t, index) = as_value(x)
# Constants
@inline next_level_local_ad(x::ConstantWrapper, t, index) = x
# Nested states
@inline next_level_local_ad(x::NamedTuple, E, index) = local_ad(x, index, E)

@inline function Base.getproperty(state::LocalStateAD{T, I, E}, f::Symbol) where {T, I, E}
    index = getfield(state, :index)
    inner_state = getfield(state, :data)
    val = getproperty(inner_state, f)
    return next_level_local_ad(val, E, index)
end

@inline Base.getindex(state::LocalStateAD, s::Symbol) = Base.getproperty(state, s)

@inline function Base.getproperty(state::ValueStateAD{T}, f::Symbol) where {T}
    inner_state = getfield(state, :data)
    val = getproperty(inner_state, f)
    return as_value(val)
end

"""
    local_ad(state::T, index::I, ad_tag::∂T) where {T, I<:Integer, ∂T}

Create local_ad for state for index I of AD tag of type ad_tag
"""
@inline local_ad(state, index, ad_tag) = local_state_ad(state, index, ad_tag)

@inline function local_state_ad(state::T, index::I, ad_tag::∂T) where {T, I<:Integer, ∂T}
    return LocalStateAD{T, I, ad_tag}(index, state)
end

function Base.show(io::IO, t::MIME"text/plain", x::LocalStateAD{T, I, E}) where {T, I, E}
    print(io, "Local state for $(unpack_tag(E)) -> $(getfield(x, :index)) with fields $(keys(getfield(x, :data)))")
end
