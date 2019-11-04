export DiaSymSemiseparableChol, trinv

struct DiaSymSemiseparableChol <: SymSemiseparableCholesky
    n::Int64
    p::Int64
    U::AbstractArray
    W::AbstractArray
    ds::AbstractArray
end

# Constuctors
DiaSymSemiseparableChol(U::AbstractArray, W::AbstractArray, ds::AbstractArray) = DiaSymSemiseparableChol(size(U,1),size(U,2),U,W,ds);

# Mappings
mul!(y::AbstractArray, L::DiaSymSemiseparableChol, 		            b::AbstractArray) =   dss_tri_mul!(y, L.U,   L.W,   L.ds,   b);
mul!(y::AbstractArray, L::AdjointOperator{DiaSymSemiseparableChol}, b::AbstractArray) =  dssa_tri_mul!(y, L.A.U, L.A.W, L.A.ds, b);
inv!(y::AbstractArray, L::DiaSymSemiseparableChol, 		      		b::AbstractArray) =   dss_forward!(y, L.U,   L.W,   L.ds,   b);
inv!(y::AbstractArray, L::AdjointOperator{DiaSymSemiseparableChol}, b::AbstractArray) =  dssa_backward!(y, L.A.U, L.A.W, L.A.ds, b);
newlogdet(L::DiaSymSemiseparableChol) = dss_logdet(L.ds)
newlogdet(L::AdjointOperator{DiaSymSemiseparableChol}) = dss_logdet(L.A.ds)

#### Inverse of a EGRQSCholesky using ####
function inv(L::DiaSymSemiseparableChol)
	return L'\(L\Diagonal(ones(L.n)))
end
function inv(L::DiaSymSemiseparableChol, b::AbstractArray)
	return L'\(L\b)
end

# Traces, norms and determinants
function fro_norm_L(L::DiaSymSemiseparableChol)
	return sum(squared_norm_cols(L.U, L.W, L.ds))
end

function trinv(L::DiaSymSemiseparableChol)
	dbar = L.ds;
	Y, Z = dss_create_yz(L.U, L.W, dbar)
	return sum(squared_norm_cols(Y,Z, dbar.^(-1)))
end
function tr(L::DiaSymSemiseparableChol)
	return sum(L.ds)
end
function tr(Ky::DiaSymSemiseparableChol, K::SymSemiseparable)
	n = Ky.n;
	p = Ky.p;
	c = Ky.ds;
	U = K.U;
	V = K.V;
	Y, Z = dss_create_yz(Ky.U, Ky.W, Ky.ds);
	b = 0;
	P = zeros(p,p);
	R = zeros(p,p);
	for k = 1:Ky.n
		yk = Y[k,:];
		zk = Z[k,:];
		uk = U[k,:];
		vk = V[k,:];
		cki = c[k]^(-1);
		b += yk'*P*yk + 2*yk'*R*uk*cki + uk'*vk*(cki^2);
		P += ((uk'*vk)*zk)*zk' + zk*(R*uk)' + (R*uk)*zk';
		R += zk*vk';
	end
	return b
end


########################################################################
#### Cholesky factorization of Higher-order quasiseparable matrices ####
########################################################################

#### Creating W and ds s.t. L = tril(UW',-1) + diag(ds) ####
function dss_create_wdbar(U::AbstractArray, V::AbstractArray, d::AbstractArray)
    n, m = size(U);
    P  = zeros(m, m);
    W  = zeros(n, m);
    dbar = zeros(n);
    for i = 1:n
        tmpU  = U[i,:]
        tmpW  = V[i,:] - P*tmpU;
        tmpds = sqrt(tmpU'*tmpW + d[i]);
        tmpW  = tmpW/tmpds
        W[i,:] = tmpW;
        dbar[i] = tmpds;
        P += tmpW*tmpW';
    end
    return W, dbar
end
#### Multiplying with Ld ####
function dss_tri_mul!(Y::AbstractArray,U::AbstractArray,W::AbstractArray,
                     ds::AbstractArray,X::AbstractArray)
     n, m = size(U)
     mx = size(X, 2)
     Wbar = zeros(m, mx)
     for i = 1:n
         tmpW = W[i,:]
         tmpU = U[i,:]
         tmpX = X[i:i,:];
         Y[i,:] = tmpU'*Wbar + ds[i]*tmpX;
         Wbar  +=  tmpW*tmpX;
     end
end
#### Adjoint of Ld ####
function dssa_tri_mul!(Y::AbstractArray,U::AbstractArray,W::AbstractArray,
                      ds::AbstractArray,X::AbstractArray)
     n, m = size(U)
     mx = size(X, 2)
     Ubar = zeros(m,mx)
     for i = n:-1:1
         tmpW = W[i,:]
         tmpU = U[i,:]
         tmpX = X[i:i,:]
         Y[i,:] = tmpW'*Ubar + ds[i]*tmpX;
         Ubar = Ubar + tmpU*tmpX;
     end
end
#### Forward substitution ####
function dss_forward!(X::AbstractArray,U::AbstractArray,W::AbstractArray,
                     ds::AbstractArray,B::AbstractArray)
    n, m = size(U)
    mx = size(B,2)
    Wbar = zeros(m,mx);
    for i = 1:n
        tmpU = U[i,:];
        tmpW = W[i,:];
        X[i:i,:] = (B[i:i,:] - tmpU'*Wbar)/ds[i];
        Wbar += tmpW .* X[i:i,:];
    end
end
#### Backward substitution ####
function dssa_backward!(X::AbstractArray,U::AbstractArray,W::AbstractArray,
                       ds::AbstractArray,B::AbstractArray)
    n, m = size(U)
    mx = size(B,2)
    Ubar = zeros(m,mx);
    for i = n:-1:1
        tmpU = U[i,:];
        tmpW = W[i,:];
        X[i:i,:] = (B[i:i,:] - tmpW'*Ubar)/ds[i];
        Ubar += tmpU .* X[i:i,:];
    end
end

#### Squared norm of columns of L = tril(UW',-1) + diag(dbar) ####
function squared_norm_cols(U::AbstractArray,W::AbstractArray,
                        dbar::AbstractArray)
    n, m = size(U)
    P = zeros(m, m)
    c = zeros(n)
    for i = n:-1:1
        tmpW = W[i,:]
        tmpU = U[i,:]
        c[i]  = dbar[i]^2 + tmpW'*P*tmpW
        P += tmpU*tmpU'
    end
    return c
end
#### Implicit inverse of  L = tril(UW',-1) + diag(dbar) ####
function dss_create_yz(U::AbstractArray, W::AbstractArray,
                    dbar::AbstractArray)
    n, m = size(U)
    Y = zeros(n,m)
    Z = zeros(n,m)
    dss_forward!(Y, U, W, dbar, U)
    dssa_backward!(Z, U, W, dbar, W)
    # Probably best not to use inv
    return Y, Z*inv(U'*Z - I)
end
#### Log-determinant ####
function dss_logdet(d::AbstractArray)
    return sum(log.(d))
end
