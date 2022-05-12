export SymSemiseparableChol

struct SymSemiseparableChol <: SymSemiseparableCholesky
    n::Int64
    p::Int64
    U::AbstractArray
    W::AbstractArray
end

# Constuctors
function SymSemiseparableChol(U::AbstractArray, W::AbstractArray)
	if size(U,1) == size(W,1) && size(U,2) == size(W,2)
		return SymSemiseparableChol(size(U,1),size(U,2),U,W)
	else
		error("Dimension mismatch between the generators U and W")
	end
end

# Mappings
mul!(y::AbstractArray, L::SymSemiseparableChol, 		         b::AbstractArray) =
	ss_tri_mul!(y, L.U,   L.W,   b)
mul!(y::AbstractArray, L::AdjointOperator{SymSemiseparableChol}, b::AbstractArray) =
	ssa_tri_mul!(y, L.A.U, L.A.W, b)
inv!(y::AbstractArray, L::SymSemiseparableChol, 		    	 b::AbstractArray) =
	ss_forward!(y, L.U,   L.W,   b)
inv!(y::AbstractArray, L::AdjointOperator{SymSemiseparableChol}, b::AbstractArray) =
	ssa_backward!(y, L.A.U, L.A.W, b)
newlogdet(L::SymSemiseparableChol) = ss_logdet(L.U, L.W)
newlogdet(L::AdjointOperator{SymSemiseparableChol}) = ss_logdet(L.A.U, L.A.W)

#### Inverse of a SymSemiseparableChol using ####
function inv(L::SymSemiseparableChol)
	return L'\(L\Diagonal(ones(L.n)))
end
function inv(L::SymSemiseparableChol, b::AbstractArray)
	return L'\(L\b)
end

########################################################################
#### Cholesky factoriaztion of:                                     ####
####    Extended generator representable {p}-semiseperable matrices ####
########################################################################

# Creating W, s.t. L = tril(UW')
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
#### Triangular product ####
function ss_tri_mul!(Y::AbstractArray,U::AbstractArray,
                     W::AbstractArray,X::AbstractArray)
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
#### Adjoint triangular product ####
function ssa_tri_mul!(Y::AbstractArray,U::AbstractArray,
                      W::AbstractArray,X::AbstractArray)
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
#### Forward substitution (solve Lx = b) ####
function ss_forward!(X::AbstractArray, U::AbstractArray,
                     W::AbstractArray, B::AbstractArray)
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
#### Backward substitution (solve L'x = b) ####
function ssa_backward!(X::AbstractArray, U::AbstractArray,
                       W::AbstractArray, B::AbstractArray)
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
function ss_logdet(U::AbstractArray, W::AbstractArray)
    n = size(U, 1)
    a = 0
    for i = 1:n
        a += log(dot(U[i,:],W[i,:]))
    end
    return a
end
