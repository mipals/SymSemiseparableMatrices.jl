module SymSemiseparableMatrices

# Importing Relevant Packages
using LinearAlgebra

# Must be imported here, as it is relevant before "syntax.jl"
import LinearAlgebra: inv!, tr, mul!, logdet, ldiv!, cholesky, det
import Base: inv, size, eltype, adjoint, *, \, size, eltype, Matrix, getindex, getproperty

# Creating Types
struct SymSemiseparableMatrix{T,UU<:AbstractMatrix{T},VV<:AbstractMatrix{T}} <: AbstractMatrix{T}
    n::Int64
    p::Int64
    Ut::UU
    Vt::VV
end
struct SymSemiseparableCholesky{T,UU<:AbstractMatrix{T},WW<:AbstractMatrix{T}} <: AbstractMatrix{T}
    n::Int64
    p::Int64
    Ut::UU
    Wt::WW
end
struct DiaSymSemiseparableMatrix{T,UU<:AbstractMatrix{T},VV<:AbstractMatrix{T},dd<:AbstractVector{T}} <: AbstractMatrix{T}
    n::Int64
    p::Int64
    Ut::UU
    Vt::VV
    d::dd
end
struct DiaSymSemiseparableCholesky{T,UU<:AbstractMatrix{T},WW<:AbstractMatrix{T},dd<:AbstractVector{T}} <: AbstractMatrix{T}
    n::Int64
    p::Int64
    Ut::UU
    Wt::WW
    d::dd
end

# Properties
# include("adjointoperator.jl")

# Matrices
include("SymSemiseparableMatrix.jl")
include("SymSemiseparableCholesky.jl")
include("DiaSymSemiseparableMatrix.jl")
include("DiaSymSemiseparableCholesky.jl")

# Operator overloadings
# include("operator_overloading.jl")

# Spline Kernel
include("spline_kernel.jl")

# Exporting Relevant 
export SymSemiseparableMatrix, SymSemiseparableCholesky
export DiaSymSemiseparableMatrix, DiaSymSemiseparableCholesky
export trinv

end # module
