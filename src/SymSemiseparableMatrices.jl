module SymSemiseparableMatrices
using LinearAlgebra

abstract type SymSemiseparableMatrix end
abstract type SymSemiseparableCholesky  <: SymSemiseparableMatrix end

# Must be imported here, as it is relevant before "syntax.jl"
import LinearAlgebra: inv!, tr, mul!
import Base: inv, size, eltype

export SymSemiseparableMatrix,
       SymSemiseparableCholesky



# Properties
include("adjointoperator.jl")

# Matrices
include("matrices/symsep.jl")
include("matrices/symsepchol.jl")
include("matrices/diasymsep.jl")
include("matrices/diasymsepchol.jl")

include("syntax.jl")

# More constructors
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

# # 2D tensor algorithms
# include("algorithms/tensor_kernel.jl")
# include("operators/tensorkernel.jl")
# include("operators/tensoreigenkernel.jl")


end # module
