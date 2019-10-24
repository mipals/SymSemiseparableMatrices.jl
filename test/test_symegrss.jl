# Creating a test problem
n = 100;
p = 5;
U = randn(n,p);
V = randn(n,p);
K = SymEGRSSMatrix(U,V);
x = randn(n);
Kfull = tril(U*V') + triu(V*U',1);

# Testing multiplication
@test isapprox(K*x, Kfull*x, atol = 1e-6)
#  Testing multiplication with the adjoint operator
@test isapprox(K'*x, Kfull'*x, atol = 1e-6)
