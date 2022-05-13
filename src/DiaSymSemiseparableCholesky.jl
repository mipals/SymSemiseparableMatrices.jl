#==========================================================================================
                                Constructors
==========================================================================================#
function DiaSymSemiseparableCholesky(U::AbstractArray, W::AbstractArray, ds::AbstractArray) 
    return DiaSymSemiseparableCholesky(size(U,1),size(U,2),U,W,ds)
end
function DiaSymSemiseparableCholesky(U::AbstractArray, V::AbstractArray, σn, σf)
    n, p = size(U)
    W, dbar = dss_create_wdbar(σf*U, σf*V, ones(n)*σn^2)
    return DiaSymSemiseparableCholesky(n, p, σf*U, W, dbar)
end
function DiaSymSemiseparableCholesky(L::DiaSymSemiseparableMatrix)
    W, dbar = dss_create_wdbar(L.Ut, L.Vt, L.d)
    return DiaSymSemiseparableCholesky(L.n, L.p, L.Ut, W, dbar)
end

#==========================================================================================
                        Defining Matrix Properties
==========================================================================================#
Matrix(K::DiaSymSemiseparableCholesky) = getproperty(K,:L)
size(K::DiaSymSemiseparableCholesky) = (K.n, K.n)
size(K::DiaSymSemiseparableCholesky,d::Int) = (1 <= d && d <=2) ? size(K)[d] : throw(ArgumentError("Invalid dimension $d"))

function getindex(K::DiaSymSemiseparableCholesky, i::Int, j::Int)
	i > j && return dot(K.Ut[:,i], K.Wt[:,j])
	i == j && return K.d[i]
	return 0
end

function getproperty(K::DiaSymSemiseparableCholesky, d::Symbol)
    if d === :U
        return UpperTriangular(triu(K.Wt'*K.Ut,1) + Diagonal(K.d))
    elseif d === :L
        return LowerTriangular(tril(K.Ut'*K.Wt,-1) + Diagonal(K.d))
    else
        return getfield(K, d)
    end
end

Base.propertynames(F::DiaSymSemiseparableCholesky, private::Bool=false) =
    (:U, :L, (private ? fieldnames(typeof(F)) : ())...)

function Base.show(io::IO, mime::MIME{Symbol("text/plain")},
		 K::DiaSymSemiseparableCholesky{<:Any,<:AbstractArray,<:AbstractArray,<:AbstractArray})
    summary(io, K); println(io)
    show(io, mime, K.L)
end
#==========================================================================================
                    Defining multiplication and inverse
==========================================================================================#
function mul!(y::AbstractArray, L::DiaSymSemiseparableCholesky, b::AbstractArray) 
    dss_tri_mul!(y, L.Ut, L.Wt, L.d, b)
    return y
end
function mul!(y::AbstractArray, L::Adjoint{<:Any,<:DiaSymSemiseparableCholesky}, b::AbstractArray)
    dssa_tri_mul!(y, L.parent.Ut, L.parent.Wt, L.parent.d, b)
    return y
end
function inv!(y::AbstractArray, L::DiaSymSemiseparableCholesky, b::AbstractArray) 
    dss_forward!(y, L.Ut, L.Wt, L.d, b)
    return y
end
function inv!(y::AbstractArray, L::Adjoint{<:Any,<:DiaSymSemiseparableCholesky}, b::AbstractArray) 
    dssa_backward!(y, L.parent.Ut, L.parent.Wt, L.parent.d, b)
    return y
end

#### Inverse of a EGRQSCholesky using the Cholesky factorization ####
function inv(L::DiaSymSemiseparableCholesky)
	return L'\(L\Diagonal(ones(L.n)))
end
function inv(L::DiaSymSemiseparableCholesky, b::AbstractArray)
	return L'\(L\b)
end

#==========================================================================================
                        Relevant Linear Algebra Routines
==========================================================================================#
function fro_norm_L(L::DiaSymSemiseparableCholesky)
	return sum(squared_norm_cols(L.Ut, L.Wt, L.d))
end

function trinv(L::DiaSymSemiseparableCholesky)
	dbar = L.d
	Y, Z = dss_create_yz(L.Ut, L.Wt, dbar)
	return sum(squared_norm_cols(Y,Z, dbar.^(-1)))
end
function tr(L::DiaSymSemiseparableCholesky)
	return sum(L.d)
end
function tr(Ky::DiaSymSemiseparableCholesky, K::SymSemiseparableMatrix)
	p = Ky.p
	c = Ky.d
	U = K.Ut
	V = K.Vt
	Y, Z = dss_create_yz(Ky.Ut, Ky.Wt, Ky.d)
	b = 0.0
	P = zeros(p,p)
	R = zeros(p,p)
	@inbounds for k = 1:Ky.n
		yk = @view Y[k,:]
		zk = @view Z[k,:]
		uk = @view U[k,:]
		vk = @view V[k,:]
		cki = c[k]^(-1)
		b += yk'*P*yk + 2*yk'*R*uk*cki + uk'*vk*(cki^2)
		P += ((uk'*vk)*zk)*zk' + zk*(R*uk)' + (R*uk)*zk'
		R += zk*vk';
	end
	return b
end
#===========================================================================================
                Choleskyesky factorization of Higher-order quasiseparable matrices
===========================================================================================#
"""
    dss_create_wdbar(U, V, d)

Computes `W` and `dbar` such that, `L = tril(UW',-1) + diag(ds)`.
"""
function dss_create_wdbar(U, V, d)
    m, n = size(U)
    P  = zeros(m, m)
    W  = zeros(m, n)
    dbar = zeros(n)
    for i = 1:n
        tmpU  = @view U[:,i]
        tmpW  = V[:,i] - P*tmpU
        tmpds = sqrt(tmpU'*tmpW + d[i])
        tmpW  = tmpW/tmpds
        W[:,i] = tmpW
        dbar[i] = tmpds
        P += tmpW*tmpW'
    end
    return W, dbar
end
#### Multiplying with Ld ####
function dss_tri_mul!(Y,U,W,ds,X)
    m, n = size(U)
    mx = size(X, 2)
    Wbar = zeros(m, mx)
    @inbounds for i = 1:n
        tmpW = @view W[:,i]
        tmpU = @view U[:,i]
        tmpX = @view X[i:i,:]
        Y[i,:] = tmpU'*Wbar + ds[i]*tmpX
        Wbar  += tmpW*tmpX
    end
end
#### Adjoint of Ld ####
function dssa_tri_mul!(Y,U,W,ds,X)
    m, n = size(U)
    mx = size(X, 2)
    Ubar = zeros(m,mx)
    @inbounds for i = n:-1:1
        tmpW = @view W[:,i]
        tmpU = @view U[:,i]
        tmpX = @view X[i:i,:]
        Y[i,:] = tmpW'*Ubar + ds[i]*tmpX;
        Ubar = Ubar + tmpU*tmpX;
    end
end
#### Forward substitution ####
function dss_forward!(X,U,W,ds,B)
    m, n = size(U)
    mx = size(B,2)
    Wbar = zeros(m,mx)
    @inbounds for i = 1:n
        tmpU = @view U[:,i]
        tmpW = @view W[:,i]
        X[i:i,:] = (B[i:i,:] - tmpU'*Wbar)/ds[i]
        Wbar += tmpW .* X[i:i,:]
    end
end
#### Backward substitution ####
function dssa_backward!(X,U,W,ds,B)
    m, n = size(U)
    mx = size(B,2)
    Ubar = zeros(m,mx)
    @inbounds for i = n:-1:1
        tmpU = @view U[:,i]
        tmpW = @view W[:,i]
        X[i:i,:] = (B[i:i,:] - tmpW'*Ubar)/ds[i]
        Ubar += tmpU .* X[i:i,:]
    end
end

#### Squared norm of columns of L = tril(UW',-1) + diag(dbar) ####
function squared_norm_cols(U,W,dbar)
    m, n = size(U)
    P = zeros(m, m)
    c = zeros(n)
    @inbounds for i = n:-1:1
        tmpW = @view W[:,i]
        tmpU = @view U[:,i]
        c[i]  = dbar[i]^2 + tmpW'*P*tmpW
        P += tmpU*tmpU'
    end
    return c
end
#### Implicit inverse of  L = tril(UW',-1) + diag(dbar) ####
function dss_create_yz(U, W,dbar)
    m, n = size(U)
    Y = zeros(n,m)
    Z = zeros(n,m)
    dss_forward!(Y, U, W, dbar, U)
    dssa_backward!(Z, U, W, dbar, W)
    # Probably best not to use inv
    return Y, Z*inv(U'*Z - I)
end
#### Log-determinant ####
dss_logdet(d) = sum(log,d)
newlogdet(L::DiaSymSemiseparableCholesky) = dss_logdet(L.d)
newlogdet(L::Adjoint{<:Any,<:DiaSymSemiseparableCholesky}) = dss_logdet(L.parent.d)
