#==========================================================================================
                                Constructors
==========================================================================================#
function SymSemiseparable(U::AbstractArray, V::AbstractArray)
    if size(U,1) == size(V,1) && size(U,2) == size(V,2)
        return SymSemiseparable(size(U,1),size(U,2),U,V);
    else
        error("Dimension mismatch between generators U and V")
    end
end
function SymSemiseparable(L::SymSemiseparableCholesky) 
    return SymSemiseparable(L.n, L.p, L.U, ss_create_v(L.U, L.W))
end
#==========================================================================================
                        Defining multiplication and inverse
==========================================================================================#
mul!(y::AbstractArray, L::SymSemiseparable, x::AbstractArray) =
    ss_mul_mat!(y, L.U, L.V, x)
mul!(y::AbstractArray, L::AdjointOperator{SymSemiseparable}, x::AbstractArray) =
    ss_mul_mat!(y, L.A.U, L.A.V, x)
function inv!(y, K::SymSemiseparable, b::AbstractArray)
	L = SymSemiseparableCholesky(K)
	y[:,:] = L'\(L\b)
end
function inv!(y, K::AdjointOperator{SymSemiseparable}, b::AbstractArray)
	L = SymSemiseparableCholesky(K.A)
	y[:,:] = L'\(L\b)
end
#==========================================================================================
                        Defining the Linear Algebra
==========================================================================================#
"""
    ss_mul_mat!(Y, U, V, X)

Computes the matrix-matrix product `Y = (tril(U*V') + triu(V*U',1))*X` in linear complexity.
"""
function ss_mul_mat!(Y, U, V, X)
    n, m = size(U)
    mx = size(X,2)
    Vbar = zeros(m,mx)
    Ubar = U'*X
    @inbounds for i = 1:n
        tmpV = @view V[i,:]
        tmpU = @view U[i,:]
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
    n,m = size(U)
    V = zeros(n,m)
    P = zeros(m,m)
    @inbounds for i = 1:n
        tmpW = @view W[i,:]
        tmpU = @view U[i,:]
        tmpV = tmpW*(tmpU'*tmpW)
        V[i,:] = tmpV + P*tmpU
        P += tmpW*tmpW'
    end
    return V
end
