#==========================================================================================
                                Struct & Constructors
==========================================================================================#
function DiaSymSemiseparableCholesky(U::AbstractArray, W::AbstractArray, ds::AbstractArray) 
    return DiaSymSemiseparableCholesky(size(U,1),size(U,2),U,W,ds)
end
function DiaSymSemiseparableCholesky(U::AbstractArray, V::AbstractArray, σn, σf)
    n, p = size(U)
    W, dbar = dss_create_wdbar(σf*U, σf*V, ones(n)*σn^2)
    return DiaSymSemiseparableCholesky(n, p, σf*U, W, dbar)
end
function DiaSymSemiseparableCholesky(L::DiaSymSemiseparable)
    W, dbar = dss_create_wdbar(L.U, L.V, L.d)
    return DiaSymSemiseparableCholesky(L.n, L.p, L.U, W, dbar)
end
#==========================================================================================
                    Defining multiplication and inverse
==========================================================================================#
function mul!(y::AbstractArray, L::DiaSymSemiseparableCholesky, b::AbstractArray) 
    dss_tri_mul!(y, L.U, L.W, L.ds, b)
end
function mul!(y::AbstractArray, L::AdjointOperator{DiaSymSemiseparableCholesky}, b::AbstractArray)
    dssa_tri_mul!(y, L.A.U, L.A.W, L.A.ds, b)
end
function inv!(y::AbstractArray, L::DiaSymSemiseparableCholesky, b::AbstractArray) 
    dss_forward!(y, L.U, L.W, L.ds, b)
end
function inv!(y::AbstractArray, L::AdjointOperator{DiaSymSemiseparableCholesky}, b::AbstractArray) 
    dssa_backward!(y, L.A.U, L.A.W, L.A.ds, b)
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
	return sum(squared_norm_cols(L.U, L.W, L.ds))
end

function trinv(L::DiaSymSemiseparableCholesky)
	dbar = L.ds;
	Y, Z = dss_create_yz(L.U, L.W, dbar)
	return sum(squared_norm_cols(Y,Z, dbar.^(-1)))
end
function tr(L::DiaSymSemiseparableCholesky)
	return sum(L.ds)
end
function tr(Ky::DiaSymSemiseparableCholesky, K::SymSemiseparable)
	n = Ky.n
	p = Ky.p
	c = Ky.ds
	U = K.U
	V = K.V
	Y, Z = dss_create_yz(Ky.U, Ky.W, Ky.ds)
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
newlogdet(L::DiaSymSemiseparableCholesky) = dss_logdet(L.ds)
newlogdet(L::AdjointOperator{DiaSymSemiseparableCholesky}) = dss_logdet(L.A.ds)
#===========================================================================================
                Choleskyesky factorization of Higher-order quasiseparable matrices
===========================================================================================#
#### Creating W and ds s.t. L = tril(UW',-1) + diag(ds) ####
"""
    dss_create_wdbar(U, V, d)

Computes `W` and `dbar` such that, `L = tril(UW',-1) + diag(ds)`.
"""
function dss_create_wdbar(U, V, d)
    n, m = size(U)
    P  = zeros(m, m)
    W  = zeros(n, m)
    dbar = zeros(n)
    for i = 1:n
        tmpU  = @view U[i,:]
        tmpW  = V[i,:] - P*tmpU
        tmpds = sqrt(tmpU'*tmpW + d[i])
        tmpW  = tmpW/tmpds
        W[i,:] = tmpW
        dbar[i] = tmpds
        P += tmpW*tmpW'
    end
    return W, dbar
end
#### Multiplying with Ld ####
function dss_tri_mul!(Y,U,W,ds,X)
    n, m = size(U)
    mx = size(X, 2)
    Wbar = zeros(m, mx)
    @inbounds for i = 1:n
        tmpW = @view W[i,:]
        tmpU = @view U[i,:]
        tmpX = @view X[i:i,:]
        Y[i,:] = tmpU'*Wbar + ds[i]*tmpX;
        Wbar  +=  tmpW*tmpX;
    end
end
#### Adjoint of Ld ####
function dssa_tri_mul!(Y,U,W,ds,X)
    n, m = size(U)
    mx = size(X, 2)
    Ubar = zeros(m,mx)
    @inbounds for i = n:-1:1
        tmpW = @view W[i,:]
        tmpU = @view U[i,:]
        tmpX = @view X[i:i,:]
        Y[i,:] = tmpW'*Ubar + ds[i]*tmpX;
        Ubar = Ubar + tmpU*tmpX;
    end
end
#### Forward substitution ####
function dss_forward!(X,U,W,ds,B)
    n, m = size(U)
    mx = size(B,2)
    Wbar = zeros(m,mx)
    @inbounds for i = 1:n
        tmpU = @view U[i,:]
        tmpW = @view W[i,:]
        X[i:i,:] = (B[i:i,:] - tmpU'*Wbar)/ds[i]
        Wbar += tmpW .* X[i:i,:]
    end
end
#### Backward substitution ####
function dssa_backward!(X,U,W,ds,B)
    n, m = size(U)
    mx = size(B,2)
    Ubar = zeros(m,mx)
    @inbounds for i = n:-1:1
        tmpU = @view U[i,:]
        tmpW = @view W[i,:]
        X[i:i,:] = (B[i:i,:] - tmpW'*Ubar)/ds[i]
        Ubar += tmpU .* X[i:i,:]
    end
end

#### Squared norm of columns of L = tril(UW',-1) + diag(dbar) ####
function squared_norm_cols(U,W,dbar)
    n, m = size(U)
    P = zeros(m, m)
    c = zeros(n)
    @inbounds for i = n:-1:1
        tmpW = @view W[i,:]
        tmpU = @view U[i,:]
        c[i]  = dbar[i]^2 + tmpW'*P*tmpW
        P += tmpU*tmpU'
    end
    return c
end
#### Implicit inverse of  L = tril(UW',-1) + diag(dbar) ####
function dss_create_yz(U, W,dbar)
    n, m = size(U)
    Y = zeros(n,m)
    Z = zeros(n,m)
    dss_forward!(Y, U, W, dbar, U)
    dssa_backward!(Z, U, W, dbar, W)
    # Probably best not to use inv
    return Y, Z*inv(U'*Z - I)
end
#### Log-determinant ####
function dss_logdet(d)
    return sum(log.(d))
end
