using SymSemiseparableMatrices
using LinearAlgebra
include("spline_kernel.jl")

U, V = spline_kernel(Vector(0.1:0.01:1), 2); # Creating Input Generators resulting in a PSD K
K = SymSemiseparable(U,V); # Generator symmetric semiseparable matrix
x = randn(K.n); # Test vector

K*x

K'*x

C = SymSemiseparableChol(K); # Creating Cholesky factorization
C*x
C'*x




function ssa_tri_mul!(Y::AbstractArray,U::AbstractArray,
                      W::AbstractArray,X::AbstractArray)
     n, m = size(U);
     mx = size(X,2);
     Ubar = zeros(m,mx);
     Ubar = Ubar + U'*X;
     for i = 1:n
         tmpW = W[i,:]
         tmpU = U[i,:]
         tmpX = X[i:i,:]
         Y[i,:] = Ubar'*tmpW;
         Ubar  -= tmpU .* tmpX
     end
end
X = x
Y = zeros(size(x));
W = C.W
ssa_tri_mul!(Y,U,W,X)


function ss_tri_mul!(Y::AbstractArray,U::AbstractArray,
                     W::AbstractArray,X::AbstractArray)
     n, m = size(U);
     mx = size(X,2);
     Wbar = zeros(m,mx);
     for i = 1:n
         tmpW   = W[i,:]
         tmpU   = U[i,:]
         tmpX   = X[i:i,:]
         #Vbar += tmpV .* X[i:i,:];
         Wbar  += tmpW .* tmpX;
         Y[i,:] = Wbar'*tmpU;
     end
end
ss_tri_mul!(Y,U,W,X)
