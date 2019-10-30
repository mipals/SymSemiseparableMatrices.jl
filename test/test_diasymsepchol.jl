include("spline_kernel.jl")

# Removing t = 0, such that Σ is invertible
t = Vector(0.1:0.1:1)
n = length(t);

# Creating a test matrix Σ = tril(UV') + triu(VU',1) that is PSD
p = 2;
U, V = spline_kernel(t, p);
Σ    = spline_kernel_matrix(U, V) + I;
chol = cholesky(Σ)

# Creating a symmetric exended generator representable semiseperable matrix
K  = DiaSymSemiseparable(U,V,ones(n))
# Calculating its Cholesky factorization
L = DiaSymSemiseparableChol(K)
# Creating a test vector
x = randn(n);

# Testing multiplication
@test isapprox(L*x, chol.L*x, atol=1e-6)
@test isapprox(L'*x, chol.U*x, atol=1e-6)

# Testing inverses (Using Cholesky factorizations)
@test isapprox(L'\(L\x), chol.U\(chol.L\x), atol=1e-6)
@test isapprox(inv(L), chol.U\(chol.L\Diagonal(ones(K.n))), atol=1e-6)
@test isapprox(inv(L,x), chol.U\(chol.L\x), atol=1e-6)
@test isapprox(K*(K\x), x)
@test isapprox(K'*(K'\x), x)

# Testing logdet
@test isapprox(logdet(L), logdet(chol.L), atol=1e-10)
@test isapprox(logdet(L'), logdet(chol.U), atol=1e-10)

# Testing trace
@test isapprox(tr(L), tr(chol.L))




#@test isapprox(fro_norm_L(L),tr(Σ))


# Testing going back
#Kl = DiaSymSemiseparable(L)
#@test isapprox(Kl*x, K*x)
