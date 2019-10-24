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
include("matrices/symegrss.jl")
# include("matricess/SemisepChol.jl")
# include("matrices/DiagSemisep.jl")
# include("matrices/DiagSemisepChol.jl")

include("syntax.jl")

# # More constructors
# EGRSSCholesky(K::EGRSSMatrix) = EGRSSCholesky(K.n, K.p, K.U, ss_create_w(K.U, K.V))
# EGRSSMatrix(L::EGRSSCholesky) = EGRSSMatrix(  L.n, L.p, L.U, ss_create_v(L.U, L.W))
# EGRQSMatrix(L::EGRSSMatrix, d::AbstractArray) = EGRQSMatrix(L.n, L.p, L.U, L.V, d)
# function EGRQSCholesky(L::EGRQSMatrix)
#       W, dbar = dss_create_wdbar(L.U, L.V, L.d)
#       EGRQSCholesky(L.n, L.p, L.U, W, dbar)
# end
# function EGRQSMatrix(L::EGRQSCholesky)
#       V, d = dss_create_vd(L.U, L.W, L.ds);
#       EGRQSMatrix(  L.n, L.p, L.U, V, d)
# end
# function EGRQSCholesky(U::AbstractArray, V::AbstractArray, σn, σf)
#       n, p = size(U)
#       W, dbar = dss_create_wdbar(σf*U, σf*V, ones(n)*σn^2)
#       EGRQSCholesky(n, p, σf*U, W, dbar)
# end
#
# # Syntax
#
# # 2D tensor algorithms
# include("algorithms/tensor_kernel.jl")
# include("operators/tensorkernel.jl")
# include("operators/tensoreigenkernel.jl")


end # module
