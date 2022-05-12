#==========================================================================================
                                Struct & Constructors
==========================================================================================#
function DiaSymSemiseparable(U::AbstractArray, V::AbstractArray, d::AbstractArray)
    if size(U,1) == size(V,1) && size(U,2) == size(V,2) && length(d) == size(U,1)
        return DiaSymSemiseparable(size(U,1),size(U,2),U,V,d);
    else
        error("Dimension mismatch between the generators U, V and d")
    end
end
function DiaSymSemiseparable(L::SymSemiseparable, d::AbstractArray)
    return DiaSymSemiseparable(L.n, L.p, L.U, L.V, d)
end
function DiaSymSemiseparable(L::DiaSymSemiseparableCholesky)
    V, d = dss_create_vd(L.U, L.W, L.ds);
    return DiaSymSemiseparable(L.n, L.p, L.U, V, d)
end

#==========================================================================================
                        Defining multiplication and inverse
==========================================================================================#
function mul!(y::AbstractArray, L::DiaSymSemiseparable, b::AbstractArray)
    dss_mul_mat!(y, L.U, L.V, L.d, b)
end
function mul!(y::AbstractArray, L::AdjointOperator{DiaSymSemiseparable}, b::AbstractArray)
    dss_mul_mat!(y, L.A.U, L.A.V, L.A.d, b)
end
function inv!(y, K::DiaSymSemiseparable, b::AbstractArray)
	L = DiaSymSemiseparableCholesky(K)
	y[:,:] = L'\(L\b)
end
function inv!(y, K::AdjointOperator{DiaSymSemiseparable}, b::AbstractArray)
	L = DiaSymSemiseparableCholesky(K.A)
	y[:,:] = L'\(L\b)
end

#===========================================================================================
                Cholesky factoriaztion of:  Higher-order quasiseparable matrices
===========================================================================================#
"""
    dss_mul_mat!(Y, U, V, d, X)

Computes the matrix-matrix product `Y = (tril(U*V') + triu(V*U',1) + diag(d))*X`.
"""
function dss_mul_mat!(Y, U, V, d, X)
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

"""
    dss_create_vd(U, W, dbar)

Computes `V` and `d` from the Cholesky factorization.
"""
function dss_create_vd(U, W, dbar)
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
