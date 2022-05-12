module SymSemiseparableMatrices

# Importing Relevant Packages
using LinearAlgebra

# Must be imported here, as it is relevant before "syntax.jl"
import LinearAlgebra: inv!, tr, mul!, logdet
import Base: inv, size, eltype, adjoint, *, \

# Creating abstract types
abstract type SymSemiseparableMatrix end
struct SymSemiseparable <: SymSemiseparableMatrix
    n::Int64
    p::Int64
    U::AbstractArray
    V::AbstractArray
end
struct SymSemiseparableCholesky <: SymSemiseparableMatrix
    n::Int64
    p::Int64
    U::AbstractArray
    W::AbstractArray
end
struct DiaSymSemiseparable <: SymSemiseparableMatrix
    n::Int64
    p::Int64
    U::AbstractArray
    V::AbstractArray
    d::AbstractArray
end
struct DiaSymSemiseparableCholesky <: SymSemiseparableMatrix
    n::Int64
    p::Int64
    U::AbstractArray
    W::AbstractArray
    ds::AbstractArray
end

# Properties
include("adjointoperator.jl")

# Matrices
include("SymSemiseparable.jl")
include("SymSemiseparableCholesky.jl")
include("DiaSymSemiseparable.jl")
include("DiaSymSemiseparableCholesky.jl")

# Operator overloadings
include("operator_overloading.jl")

# Spline Kernel
include("spline_kernel.jl")

# Exporting Relevant 
export SymSemiseparableMatrix
export SymSemiseparable, SymSemiseparableCholesky
export DiaSymSemiseparable, DiaSymSemiseparableCholesky
export trinv

end # module
