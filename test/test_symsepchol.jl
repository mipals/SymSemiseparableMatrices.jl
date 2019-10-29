include("spline_kernel.jl")

# Removing t = 0, such that Σ is invertible
t = Vector(0.1:0.1:1)
n = length(t); p = 2;

# Creating generators U,V that result in a positive-definite matrix Σ
U, V = spline_kernel(t, p);
Σ    = spline_kernel_matrix(U, V);
chol = cholesky(Σ)

# Creating a symmetric extended generator representable semiseperable matrix
K  = SymSemiseparable(U,V)
# Calculating its Cholesky factorization
Km = SymSemiseparableChol(K)
# Creating a test vector
xt = randn(n);

# Testing multiplication
@test isapprox(Km*xt, chol.L*xt, atol=1e-6)
@test isapprox(Km'*xt, chol.U*xt, atol=1e-6)

# Testing inverses (Using Cholesky factorizations)
@test isapprox(Km'\(Km\xt), chol.U\(chol.L\xt), atol=1e-6)
@test isapprox(inv(Km), chol.U\(chol.L\Diagonal(ones(K.n))), atol=1e-6)
@test isapprox(inv(Km,xt),  chol.U\(chol.L\xt), atol=1e-6)

# Testing logdet
@test isapprox(logdet(Km), logdet(chol.L), atol=1e-10)
@test isapprox(logdet(Km'), logdet(chol.U), atol=1e-10)
