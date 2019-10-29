include("spline_kernel.jl")

# Removing t = 0, such that Σ is invertible
t = Vector(0.1:0.1:1); p = 2;

# Creating generators U,V that result in a positive-definite matrix Σ
U, V = spline_kernel(t, p);
Σ    = spline_kernel_matrix(U, V);
chol = cholesky(Σ)

# Creating a symmetric extended generator representable semiseperable matrix
K = SymSemiseparable(U,V)
# Calculating its Cholesky factorization
L = SymSemiseparableChol(K)
# Creating a test vector
xt = randn(K.n);

# Testing multiplication
@test isapprox(L*xt, chol.L*xt, atol=1e-6)
@test isapprox(L'*xt, chol.U*xt, atol=1e-6)

# Testing inverses (Using Cholesky factorizations)
@test isapprox(L'\(L\xt), chol.U\(chol.L\xt), atol=1e-6)
@test isapprox(inv(L), chol.U\(chol.L\Diagonal(ones(K.n))), atol=1e-6)
@test isapprox(inv(L,xt),  chol.U\(chol.L\xt), atol=1e-6)

# Testing logdet
@test isapprox(logdet(L), logdet(chol.L), atol=1e-10)
@test isapprox(logdet(L'), logdet(chol.U), atol=1e-10)
