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
function mul!(y::AbstractArray, L::DiaSymSemiseparable, b::AbstractArray)
    dss_mul_mat!(y, L.U, L.V, L.d, b)
end
function mul!(y::AbstractArray, L::AdjointOperator{DiaSymSemiseparable}, b::AbstractArray)
    dss_mul_mat!(y, L.A.U, L.A.V, L.A.d, b)
end
function inv!(y, K::DiaSymSemiseparable, b::AbstractArray)
	L = DiaSymSemiseparableChol(K)
	y[:,:] = L'\(L\b)
end
function inv!(y, K::AdjointOperator{DiaSymSemiseparable}, b::AbstractArray)
	L = DiaSymSemiseparableChol(K.A)
	y[:,:] = L'\(L\b)
end

#===========================================================================================
                CCholesky factoriaztion of:  Higher-order quasiseparable matrices
===========================================================================================#
#### Matrix-matrix product ####
function dss_mul_mat!(Y::Array, U::Array, V::Array, d::Array, X::Array)
    n, m = size(U)
    mx = size(X,2)
    Vbar = zeros(m,mx)
    Ubar = U'*X
    @inbounds for i = 1:n
        tmpV = @view V[i,:]
        tmpU = @view U[i,:]
        tmpX = @view X[i:i,:]
        Ubar -= tmpU .* tmpX
        Vbar += tmpV .* tmpX
        Y[i,:] = tmpU'*Vbar + tmpV'*Ubar + d[i]*tmpX
    end
end
#### Going from the Cholesky factorization back to matrix Σ ####
function dss_create_vd(U::Array, W::Array, dbar::Array)
    n,m = size(U)
    d = zeros(n)
    V = zeros(n,m)
    P = zeros(m,m)
    @inbounds for i = 1:n
        tmpU = @view U[i,:]
        tmpW = @view W[i,:]
        d[i] = dbar[i]^2 - tmpU'*tmpW
        V[i,:] = P*tmpU
        P += tmpW*tmpW'
    end
    return V, d
end
