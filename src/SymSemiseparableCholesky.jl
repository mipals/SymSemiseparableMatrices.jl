#==========================================================================================
                                Constructors
==========================================================================================#
function SymSemiseparableCholesky(U::AbstractArray, W::AbstractArray)
	if size(U,1) == size(W,1) && size(U,2) == size(W,2)
		return SymSemiseparableCholesky(size(U,2),size(U,1),U,W)
	else
		error("Dimension mismatch between the generators U and W")
	end
end
function SymSemiseparableCholesky(K::SymSemiseparableMatrix)
    return SymSemiseparableCholesky(K.n, K.p, K.Ut, ss_create_w(K.Ut, K.Vt))
end
#==========================================================================================
                        Defining Matrix Properties
==========================================================================================#
Matrix(K::SymSemiseparableCholesky) = getproperty(K,:L)
size(K::SymSemiseparableCholesky) = (K.n, K.n)
size(K::SymSemiseparableCholesky, d::Int) = (1 <= d && d <=2) ? size(K)[d] : throw(ArgumentError("Invalid dimension $d"))

function getindex(K::SymSemiseparableCholesky, i::Int, j::Int)
	i >= j && return dot(K.Ut[:,i], K.Wt[:,j])
	return 0
end

function getproperty(K::SymSemiseparableCholesky, d::Symbol)
    if d === :U
        return UpperTriangular(K.Wt'*K.Ut)
    elseif d === :L
        return LowerTriangular(K.Ut'*K.Wt)
    else
        return getfield(K, d)
    end
end

# Base.propertynames(F::SymSemiseparableCholesky, private::Bool=false) =
#     (:U, :L, (private ? fieldnames(typeof(F)) : ())...)


# function Base.show(io::IO, mime::MIME{Symbol("text/plain")},
# 		 K::SymSemiseparableCholesky{<:Any,<:AbstractMatrix,<:AbstractMatrix})
#     summary(io, K); println(io)
#     show(io, mime, K.L)
# end


#==========================================================================================
                        Defining multiplication and inverse
==========================================================================================#
function mul!(y::AbstractArray, L::SymSemiseparableCholesky, b::AbstractArray)
	ss_tri_mul!(y, L.Ut,   L.Wt,   b)
    return y
end
function mul!(y::AbstractArray, L::Adjoint{<:Any,<:SymSemiseparableCholesky}, b::AbstractArray)
	ssa_tri_mul!(y, L.parent.Ut, L.parent.Wt, b)
    return y
end
function (\)(F::SymSemiseparableCholesky, B::AbstractVecOrMat)
	X = similar(B)
	ss_forward!(X,F.Ut,F.Wt,B)
	return X
end
function (\)(F::Adjoint{<:Any,<:SymSemiseparableCholesky}, B::AbstractVecOrMat)
	Y = similar(B)
	ssa_backward!(Y,F.parent.Ut,F.parent.Wt,B)
	return Y
end

newlogdet(L::SymSemiseparableCholesky) = ss_logdet(L.Ut, L.Wt)
newlogdet(L::Adjoint{<:Any,<:SymSemiseparableCholesky}) = ss_logdet(L.parent.Ut, L.parent.Wt)


function det(L::SymSemiseparableCholesky)
    dd = one(eltype(L))
    @inbounds for i in 1:L.n
        dd *= dot(L.Ut[:,i],L.Wt[:,i])
    end
    return dd
end

function logdet(L::SymSemiseparableCholesky)
    dd = zero(eltype(L))
    @inbounds for i in 1:L.n
        dd += log(dot(L.Ut[:,i],L.Wt[:,i]))
    end
    return dd
end
#### Inverse of a SymSemiseparableCholesky using ####
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
    p,n = size(U)
    wTw = zeros(p,p)
    W = zeros(p,n)
    @inbounds for j = 1:n
        tmpu = @view U[:,j]
        tmp  = V[:,j] - wTw*tmpu
        w = tmp/sqrt(abs(tmpu'*tmp))
        W[:,j] = w
        wTw = wTw + w*w'
    end
    return W
end
"""
    ss_tri_mul!(Y,U,W,X)

Computes the multiplization `tril(U*W')*X = Y` in linear complexity.
"""
function ss_tri_mul!(Y,U,W,X)
    m, n = size(U)
    mx = size(X,2)
    Wbar = zeros(m,mx)
    @inbounds for i = 1:n
        tmpW   = @view W[:,i]
        tmpU   = @view U[:,i]
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
    m, n = size(U)
    mx = size(X,2)
    Ubar = zeros(m,mx)
    Ubar = Ubar + U*X
    @inbounds for i = 1:n
        tmpW = @view W[:,i]
        tmpU = @view U[:,i]
        tmpX = @view X[i:i,:]
        Y[i,:] = Ubar'*tmpW
        Ubar  -= tmpU .* tmpX
    end
end
"""
    ss_forward!

Solves the linear system `tril(U*W')*X = B`.
"""
function ss_forward!(X,U,W, B)
    m, n = size(U)
    mx = size(B,2)
    Wbar = zeros(m, mx)
    @inbounds for i = 1:n
        tmpU = @view U[:,i]
        tmpW = @view W[:,i]
        X[i,:] = (B[i:i,:] - tmpU'*Wbar)./(tmpU'*tmpW)
        Wbar += tmpW .* X[i:i,:]
    end
end
"""
    ssa_backward!(X, U, W, B)

Solves the linear system `triu(W*U')*X = Y`.
"""
function ssa_backward!(X, U, W, B)
    m,n = size(U)
    mx = size(B,2)
    Ubar = zeros(m,mx)
    @inbounds for i = n:-1:1
        tmpU = @view U[:,i]
        tmpW = @view W[:,i]
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
