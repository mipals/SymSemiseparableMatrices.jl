n = 100;
p = 5;
U = randn(n,p);
V = randn(n,p);
K = SymEGRSSMatrix(U,V);
x = randn(n);
Kfull = tril(U*V') + triu(V*U',1);

@test isapprox(K*x, Kfull*x, atol = 1e-6)
