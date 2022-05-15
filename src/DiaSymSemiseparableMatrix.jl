#==========================================================================================
                                Constructors
==========================================================================================#
function DiaSymSemiseparableMatrix(U::AbstractArray, V::AbstractArray, d::AbstractArray)
    if size(U,1) == size(V,1) && size(U,2) == size(V,2) && length(d) == size(U,2)
        return DiaSymSemiseparableMatrix(size(U,2),size(U,1),U,V,d);
    else
        error("Dimension mismatch between the generators U, V and d")
    end
end
function DiaSymSemiseparableMatrix(L::SymSemiseparableMatrix, d::AbstractArray)
    return DiaSymSemiseparableMatrix(L.n, L.p, L.Ut, L.Vt, d)
end
function DiaSymSemiseparableMatrix(L::DiaSymSemiseparableCholesky)
    V, d = dss_create_vd(L.Ut, L.Wt, L.ds);
    return DiaSymSemiseparableMatrix(L.n, L.p, L.Ut, V, d)
end

#==========================================================================================
                        Defining Matrix Properties
==========================================================================================#
Matrix(K::DiaSymSemiseparableMatrix) = tril(K.Ut'*K.Vt) + triu(K.Vt'*K.Ut,1) + Diagonal(K.d)
size(K::DiaSymSemiseparableMatrix)   = (K.n,K.n)
LinearAlgebra.cholesky(K::DiaSymSemiseparableMatrix) = DiaSymSemiseparableCholesky(K)
function getindex(K::DiaSymSemiseparableMatrix{T}, i::Int, j::Int) where T
	i > j && return dot(K.Ut[:,i],K.Vt[:,j])
	j == i && return dot(K.Vt[:,i],K.Ut[:,j]) + K.d[i]
	return dot(K.Vt[:,i],K.Ut[:,j])
end
Base.propertynames(F::DiaSymSemiseparableMatrix, private::Bool=false) =
    (private ? fieldnames(typeof(F)) : ())

#==========================================================================================
                        Defining multiplication and inverse
==========================================================================================#
function mul!(y::AbstractArray, L::DiaSymSemiseparableMatrix, b::AbstractArray)
    dss_mul_mat!(y, L.Ut, L.Vt, L.d, b)
    return y
end
function mul!(y::AbstractArray, L::Adjoint{<:Any,<:DiaSymSemiseparableMatrix}, b::AbstractArray)
    dss_mul_mat!(y, L.parent.Ut, L.parent.Vt, L.parent.d, b)
    return y
end
function inv!(y, K::DiaSymSemiseparableMatrix, b::AbstractArray)
	L = DiaSymSemiseparableCholesky(K)
	y[:,:] = L'\(L\b)
end
function inv!(y, K::Adjoint{<:Any,<:DiaSymSemiseparableMatrix}, b::AbstractArray)
	L = DiaSymSemiseparableCholesky(K.parent)
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
    p, n = size(U)
    mx = size(X,2)
    Vbar = zeros(p,mx)
    Ubar = U*X
    @inbounds for i = 1:n
        tmpV = @view V[:,i]
        tmpU = @view U[:,i]
        tmpX = @view X[i:i,:]
        Ubar -= tmpU .* tmpX;
        Vbar += tmpV .* tmpX;
        Y[i,:] = tmpU'*Vbar + tmpV'*Ubar + d[i]*tmpX
    end
	return Y
end

"""
    dss_create_vd(U, W, dbar)

Using `L = tril(U*W',-1) + Diagonal(dbar)`, compute `V` and `d` such that
`LL = tril(UV') + triu(V'*U,1) + Diagonal(d)`.
"""
function dss_create_vd(U, W, dbar)
    m, n = size(U)
    d = zeros(n)
    V = zeros(n,m)
    P = zeros(m,m)
    @inbounds for i = 1:n
        tmpU = @view U[:,:]
        tmpW = @view W[:,:]
        d[i] = dbar[i]^2 - tmpU'*tmpW
        V[:,i] = P*tmpU
        P += tmpW*tmpW'
    end
    return V, d
end
