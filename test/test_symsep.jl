# Creating a test problem
include("spline_kernel.jl")

# Creating generators U,V that result in a positive-definite matrix K
t = Vector(0.1:0.1:1)
n = length(t); p = 2;
U, V = spline_kernel(t, p);
K = SymSemiseparable(U,V);
x = randn(K.n);
Kfull = tril(U*V') + triu(V*U',1);

# Testing multiplication
@test isapprox(K*x, Kfull*x, atol = 1e-6)

# Testing multiplication with the adjoint operator
@test isapprox(K'*x, Kfull'*x, atol = 1e-6)

# Testing linear solves
@test isapprox(K*(K\x),x, atol=1e-6)
@test isapprox(K'*(K'\x),x, atol=1e-6)
