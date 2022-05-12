module SymSemiseparableMatrices

# Importing Relevant Packages
using LinearAlgebra

# Must be imported here, as it is relevant before "syntax.jl"
import LinearAlgebra: inv!, tr, mul!, logdet
import Base: inv, size, eltype, adjoint, *, \

# Creating abstract types
abstract type SymSemiseparableMatrix end

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
export SymSemiseparableMatrix, SymSemiseparableCholesky

# More constructors - Should be moved to other files
SymSemiseparableCholesky(K::SymSemiseparable) = SymSemiseparableCholesky(K.n, K.p, K.U, ss_create_w(K.U, K.V))
SymSemiseparable(L::SymSemiseparableCholesky) = SymSemiseparable(L.n, L.p, L.U, ss_create_v(L.U, L.W))
DiaSymSemiseparable(L::SymSemiseparable, d::AbstractArray) = DiaSymSemiseparable(L.n, L.p, L.U, L.V, d)
function DiaSymSemiseparableCholesky(L::DiaSymSemiseparable)
      W, dbar = dss_create_wdbar(L.U, L.V, L.d)
      return DiaSymSemiseparableCholesky(L.n, L.p, L.U, W, dbar)
end
function DiaSymSemiseparable(L::DiaSymSemiseparableCholesky)
      V, d = dss_create_vd(L.U, L.W, L.ds);
      return DiaSymSemiseparable(L.n, L.p, L.U, V, d)
end

function DiaSymSemiseparableCholesky(U::AbstractArray, V::AbstractArray, σn, σf)
      n, p = size(U)
      W, dbar = dss_create_wdbar(σf*U, σf*V, ones(n)*σn^2)
      DiaSymSemiseparableCholesky(n, p, σf*U, W, dbar)
end


end # module
