using SymSemiseparableMatrices
using LinearAlgebra
include("spline_kernel.jl")

U, V = spline_kernel(Vector(0.1:0.01:1), 2); # Creating Input Generators resulting in a PSD K
K = SymSemiseparable(U,V); # Generator symmetric semiseparable matrix
x = ones(K.n); # Test vector

K*x

K'*x

K*(K\x)

L = SymSemiseparableChol(K); # Creating Cholesky factorization
L*x
L'*x

L\x
L'\x

Kdiag = DiaSymSemiseparable(U,V,rand(size(U,1)))
