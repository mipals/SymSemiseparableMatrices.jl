#==========================================================================================
                                Constructors
==========================================================================================#
function SymSemiseparableMatrix(U::AbstractArray, V::AbstractArray)
    if size(U,1) == size(V,1) && size(U,2) == size(V,2)
        return SymSemiseparableMatrix(size(U,2),size(U,1),U,V)
    else
        error("Dimension mismatch between generators U and V")
    end
end
function SymSemiseparableMatrix(L::SymSemiseparableCholesky) 
    return SymSemiseparableMatrix(L.n, L.p, L.Ut, ss_create_v(L.Ut, L.Wt))
end
#==========================================================================================
                        Defining Matrix Properties
==========================================================================================#
Matrix(K::SymSemiseparableMatrix)   = tril(K.Ut'*K.Vt) + triu(K.Vt'*K.Ut,1)
Base.size(K::SymSemiseparableMatrix)     = (K.n,K.n)
LinearAlgebra.cholesky(K::SymSemiseparableMatrix) = SymSemiseparableCholesky(K)
function LinearAlgebra.cholesky(K::SymSemiseparableMatrix, sigma)
	Wt, dbar = dss_create_wdbar(K.Ut,K.Vt,ones(K.n)*sigma)
	return DiaSymSemiseparableCholesky(K.n,K.p,K.Ut,Wt,dbar,)
end
function getindex(K::SymSemiseparableMatrix{T}, i::Int, j::Int) where T
	i > j && return dot(K.Ut[:,i],K.Vt[:,j])
	return dot(K.Ut[:,j],K.Vt[:,i])
end
Base.propertynames(F::SymSemiseparableMatrix, private::Bool=false) =
    (private ? fieldnames(typeof(F)) : ())

LinearAlgebra.logdet(K::SymSemiseparableMatrix) = 2.0*logdet(cholesky(K))
LinearAlgebra.det(K::SymSemiseparableMatrix)    = det(cholesky(K))^2
#==========================================================================================
                        Defining multiplication and inverse
==========================================================================================#
function LinearAlgebra.mul!(y::AbstractArray, L::SymSemiseparableMatrix, x::AbstractArray)
    ss_mul_mat!(y, L.Ut, L.Vt, x)
    return y
end
function LinearAlgebra.mul!(y::AbstractArray, L::Adjoint{<:Any,<:SymSemiseparableMatrix}, x::AbstractArray)
    ss_mul_mat!(y, L.parent.Ut, L.parent.Vt, x)
    return y
end
function LinearAlgebra.ldiv!(y::AbstractArray,K::SymSemiseparableMatrix, x::AbstractArray)
	L = cholesky(K)
    y .= L'\(L\x)
    return y
end
function LinearAlgebra.ldiv!(y::AbstractArray,K::Adjoint{<:Any,<:SymSemiseparableMatrix}, x::AbstractArray)
	L  = cholesky(K.parent)
	y .= L'\(L\x)
    return y
end
function (Base.:\)(K::SymSemiseparableMatrix, x::AbstractVecOrMat)
    y = similar(x)
	L = cholesky(K)
    y .= L'\(L\x)
    return y
end
function (Base.:\)(K::Adjoint{<:Any,<:SymSemiseparableMatrix}, x::AbstractVecOrMat)
    y = similar(x)
	L  = cholesky(K.parent)
	y .= L'\(L\x)
    return y
end
#==========================================================================================
                        Defining the Linear Algebra
==========================================================================================#
"""
    ss_mul_mat!(Y, U, V, X)

Computes the matrix-matrix product `Y = (tril(U'*V) + triu(V'*U,1))*X` in linear complexity.
"""
function ss_mul_mat!(Y, U, V, X)
    m, n = size(U)
    mx = size(X,2)
    Vbar = zeros(m,mx)
    Ubar = U*X
    @inbounds for i = 1:n
        tmpV = @view V[:,i]
        tmpU = @view U[:,i]
        Ubar -= tmpU .* X[i:i,:]
        Vbar += tmpV .* X[i:i,:]
        Y[i,:] = Vbar'*tmpU + Ubar'*tmpV
    end
end

"""
    ss_create_v(U, W)

Using `L = tril(U*W')`, compute `V` such that `LL = tril(UV') + triu(V'*U,1)`.
"""
function ss_create_v(U, W)
    m,n = size(U)
    V = zeros(n,m)
    P = zeros(m,m)
    @inbounds for i = 1:n
        tmpW = @view W[:,i]
        tmpU = @view U[:,i]
        tmpV = tmpW*(tmpU'*tmpW)
        V[i,:] = tmpV + P*tmpU
        P += tmpW*tmpW'
    end
    return V
end
