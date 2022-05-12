module SymSemiseparableMatrices

# Importing Relevant Packages
using LinearAlgebra

# Must be imported here, as it is relevant before "syntax.jl"
import LinearAlgebra: inv!, tr, mul!, logdet
import Base: inv, size, eltype, adjoint, *, \

# Creating abstract types
abstract type SymSemiseparableMatrix                              end
abstract type SymSemiseparableCholesky  <: SymSemiseparableMatrix end

# Properties
include("adjointoperator.jl")

# Matrices
include("SymSemiseparable.jl")
include("SymSemiseparableChol.jl")
include("DiaSymSemiseparable.jl")
include("DiaSymSemiseparableChol.jl")

# Operator overloadings
include("operator_overloading.jl")

# Spline Kernel
include("spline_kernel.jl")

# Exporting Relevant 
export SymSemiseparableMatrix, SymSemiseparableCholesky

# More constructors - Should be moved to other files
SymSemiseparableChol(K::SymSemiseparable) = SymSemiseparableChol(K.n, K.p, K.U, ss_create_w(K.U, K.V))
SymSemiseparable(L::SymSemiseparableChol) = SymSemiseparable(L.n, L.p, L.U, ss_create_v(L.U, L.W))
DiaSymSemiseparable(L::SymSemiseparable, d::AbstractArray) = DiaSymSemiseparable(L.n, L.p, L.U, L.V, d)
function DiaSymSemiseparableChol(L::DiaSymSemiseparable)
      W, dbar = dss_create_wdbar(L.U, L.V, L.d)
      return DiaSymSemiseparableChol(L.n, L.p, L.U, W, dbar)
end
function DiaSymSemiseparable(L::DiaSymSemiseparableChol)
      V, d = dss_create_vd(L.U, L.W, L.ds);
      return DiaSymSemiseparable(L.n, L.p, L.U, V, d)
end

function DiaSymSemiseparableChol(U::AbstractArray, V::AbstractArray, σn, σf)
      n, p = size(U)
      W, dbar = dss_create_wdbar(σf*U, σf*V, ones(n)*σn^2)
      DiaSymSemiseparableChol(n, p, σf*U, W, dbar)
end


end # module
