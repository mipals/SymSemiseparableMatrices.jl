export DiaSymSemiseparable

struct DiaSymSemiseparable <: SymSemiseparableMatrix
    n::Int64
    p::Int64
    U::AbstractArray
    V::AbstractArray
    d::AbstractArray
end

# Constuctors
function DiaSymSemiseparable(U::AbstractArray, V::AbstractArray, d::AbstractArray)
    if size(U,1) == size(V,1) && size(U,2) == size(V,2) && length(d) == size(U,1)
        return DiaSymSemiseparable(size(U,1),size(U,2),U,V,d);
    else
        error("Dimension mismatch between the generators U, V and d")
    end
end
# Mappings
mul!(y::AbstractArray, L::DiaSymSemiseparable, 		            b::AbstractArray) = dss_mul_mat!(y, L.U, L.V, L.d, b);
mul!(y::AbstractArray, L::AdjointOperator{DiaSymSemiseparable}, b::AbstractArray) = dss_mul_mat!(y, L.A.U, L.A.V, L.A.d, b);
function inv!(y, K::DiaSymSemiseparable, b::AbstractArray)
	L = DiaSymSemiseparableChol(K)
	y[:,:] = L'\(L\b)
end
function inv!(y, K::AdjointOperator{DiaSymSemiseparable}, b::AbstractArray)
	L = DiaSymSemiseparableChol(K.A)
	y[:,:] = L'\(L\b)
end

# Traces
# function tr2(K::DiaSymSemiseparable)
# 	s = 0.0;
# 	for i = 1:K.n
# 		s += K.U[i,:]'*K.V[i,:];
# 	end
# 	return s
# end

#################################################
#### Cholesky factoriaztion of:              ####
####    Higher-order quasiseparable matrices ####
#################################################

#### Matrix-matrix product ####
function dss_mul_mat!(Y::Array, U::Array, V::Array, d::Array, X::Array)
    n, m = size(U);
    mx = size(X,2);
    Vbar = zeros(m,mx);
    Ubar = U'*X;
    for i = 1:n
        tmpV = V[i,:];
        tmpU = U[i,:];
        tmpX = X[i:i,:];
        Ubar -= tmpU .* tmpX;
        Vbar += tmpV .* tmpX;
        Y[i,:] = tmpU'*Vbar + tmpV'*Ubar + d[i]*tmpX;
    end
end
#### Going from the Cholesky factorization back to matrix Σ ####
function dss_create_vd(U::Array, W::Array, dbar::Array)
    n,m = size(U)
    d = zeros(n)
    V = zeros(n,m)
    P = zeros(m,m)
    for i = 1:n
        tmpU = U[i,:]
        tmpW = W[i,:]
        d[i] = dbar[i]^2 - tmpU'*tmpW
        V[i,:] = P*tmpU
        P += tmpW*tmpW';
    end
    return V, d
end
