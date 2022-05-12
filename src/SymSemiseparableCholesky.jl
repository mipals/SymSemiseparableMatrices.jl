#==========================================================================================
                                Struct & Constructors
==========================================================================================#
function SymSemiseparableCholesky(U::AbstractArray, W::AbstractArray)
	if size(U,1) == size(W,1) && size(U,2) == size(W,2)
		return SymSemiseparableCholesky(size(U,1),size(U,2),U,W)
	else
		error("Dimension mismatch between the generators U and W")
	end
end
function SymSemiseparableCholesky(K::SymSemiseparable)
    return SymSemiseparableCholesky(K.n, K.p, K.U, ss_create_w(K.U, K.V))
end
#==========================================================================================
                        Defining multiplication and inverse
==========================================================================================#
function mul!(y::AbstractArray, L::SymSemiseparableCholesky, b::AbstractArray)
	ss_tri_mul!(y, L.U,   L.W,   b)
end
function mul!(y::AbstractArray, L::AdjointOperator{SymSemiseparableCholesky}, b::AbstractArray)
	ssa_tri_mul!(y, L.A.U, L.A.W, b)
end
function inv!(y::AbstractArray, L::SymSemiseparableCholesky, b::AbstractArray)
	ss_forward!(y, L.U,   L.W,   b)
end
function inv!(y::AbstractArray, L::AdjointOperator{SymSemiseparableCholesky}, b::AbstractArray)
	ssa_backward!(y, L.A.U, L.A.W, b)
end
newlogdet(L::SymSemiseparableCholesky) = ss_logdet(L.U, L.W)
newlogdet(L::AdjointOperator{SymSemiseparableCholesky}) = ss_logdet(L.A.U, L.A.W)

#### Inverse of a SymSemiseparableCholesky using ####
function inv(L::SymSemiseparableCholesky)
	return L'\(L\Diagonal(ones(L.n)))
end
function inv(L::SymSemiseparableCholesky, b::AbstractArray)
	return L'\(L\b)
end

#==========================================================================================
                        Defining the computations
==========================================================================================#
"""
    ss_create_w(U,V)

Creates matrix `W` such that the Cholesky factorization of the SymSemiseparable(U,V) matrix
is `L = tril(U*W')` in linear complexity.
"""
function ss_create_w(U,V)
    n,p = size(U)
    wTw = zeros(p,p)
    W = zeros(n,p)
    @inbounds for j = 1:n
        tmpu = @view U[j,:]
        tmp  = V[j,:] - wTw*tmpu;
        w = tmp/sqrt(abs(tmpu'*tmp))
        W[j,:] = w;
        wTw = wTw + w*w';
    end
    return W
end
"""
    ss_tri_mul!(Y,U,W,X)

Computes the multiplization `tril(U*W')*X = Y` in linear complexity.
"""
function ss_tri_mul!(Y,U,W,X)
    n, m = size(U)
    mx = size(X,2)
    Wbar = zeros(m,mx)
    @inbounds for i = 1:n
        tmpW   = @view W[i,:]
        tmpU   = @view U[i,:]
        tmpX   = @view X[i:i,:]
        Wbar  += tmpW .* tmpX
        Y[i,:] = Wbar'*tmpU
    end
end
"""
    ssa_tri_mul!(Y,U,W,X)

Computes the multiplization `triu(W*U')*X = Y` in linear complexity.
"""
function ssa_tri_mul!(Y,U,W,X)
    n, m = size(U)
    mx = size(X,2)
    Ubar = zeros(m,mx)
    Ubar = Ubar + U'*X
    @inbounds for i = 1:n
        tmpW = @view W[i,:]
        tmpU = @view U[i,:]
        tmpX = @view X[i:i,:]
        Y[i,:] = Ubar'*tmpW
        Ubar  -= tmpU .* tmpX
    end
end
"""
    ss_forward!

Solves the linear system `tril(U*W')*X = B`.
"""
function ss_forward!(X, U,W, B)
    n, m = size(U)
    mx = size(B,2)
    Wbar = zeros(m, mx)
    @inbounds for i = 1:n
        tmpU = @view U[i,:]
        tmpW = @view W[i,:]
        X[i,:] = (B[i:i,:] - tmpU'*Wbar)./(tmpU'*tmpW)
        Wbar += tmpW .* X[i:i,:]
    end
end
"""
    ssa_backward!(X, U, W, B)

Solves the linear system `triu(W*U')*X = Y`.
"""
function ssa_backward!(X, U, W, B)
    n, m = size(U)
    mx = size(B,2)
    Ubar = zeros(m,mx)
    @inbounds for i = n:-1:1
        tmpU = @view U[i,:]
        tmpW = @view W[i,:]
        X[i,:] = (B[i:i,:] - tmpW'*Ubar)/(tmpU'*tmpW)
        Ubar += tmpU .* X[i:i,:]
    end
end
#### Logarithm of determinant ####
function ss_logdet(U, W)
    a = 0.0
    @inbounds for i = 1:size(U, 1)
        a += log(dot(U[i,:],W[i,:]))
    end
    return a
end
