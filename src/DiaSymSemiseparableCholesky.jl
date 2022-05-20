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
                    Defining Matrix Properties / Overloading Base
==========================================================================================#
Base.Matrix(K::DiaSymSemiseparableCholesky) = getproperty(K,:L)
Base.size(K::DiaSymSemiseparableCholesky) = (K.n, K.n)
function Base.size(K::DiaSymSemiseparableCholesky,d::Int)
    return (1 <= d && d <=2) ? size(K)[d] : throw(ArgumentError("Invalid dimension $d"))
end
function Base.getindex(K::DiaSymSemiseparableCholesky, i::Int, j::Int)
	i > j && return dot(K.Ut[:,i], K.Wt[:,j])
	i == j && return K.d[i]
	return 0
end
function Base.getproperty(K::DiaSymSemiseparableCholesky, d::Symbol)
    if d === :U
        return UpperTriangular(triu(K.Wt'*K.Ut,1) + Diagonal(K.d))
    elseif d === :L
        return LowerTriangular(tril(K.Ut'*K.Wt,-1) + Diagonal(K.d))
    else
        return getfield(K, d)
    end
end
function Base.propertynames(F::DiaSymSemiseparableCholesky, private::Bool=false)
    return (:U, :L, (private ? fieldnames(typeof(F)) : ())...)
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")},
		 K::DiaSymSemiseparableCholesky{<:Any,<:AbstractArray,<:AbstractArray,<:AbstractArray})
    summary(io, K); println(io)
    show(io, mime, K.L)
end
#==========================================================================================
                            Overloading LinearAlgebra routines
==========================================================================================#
dss_logdet(d) = sum(log,d)
newlogdet(L::DiaSymSemiseparableCholesky) = dss_logdet(L.d)
newlogdet(L::Adjoint{<:Any,<:DiaSymSemiseparableCholesky}) = dss_logdet(L.parent.d)
LinearAlgebra.logdet(L::DiaSymSemiseparableCholesky) = dss_logdet(L.d)
LinearAlgebra.logdet(L::Adjoint{<:Any,<:DiaSymSemiseparableCholesky}) = dss_logdet(L.parent.d)
#==========================================================================================
                        Defining multiplication and inverse
==========================================================================================#
function LinearAlgebra.mul!(y::AbstractArray, L::DiaSymSemiseparableCholesky, b::AbstractArray)
    dss_tri_mul!(y, L.Ut, L.Wt, L.d, b)
    return y
end
function LinearAlgebra.mul!(y::AbstractArray, L::Adjoint{<:Any,<:DiaSymSemiseparableCholesky}, b::AbstractArray)
    dssa_tri_mul!(y, L.parent.Ut, L.parent.Wt, L.parent.d, b)
    return y
end
function LinearAlgebra.ldiv!(y::AbstractArray, L::DiaSymSemiseparableCholesky, b::AbstractArray)
    dss_forward!(y, L.Ut, L.Wt, L.d, b)
    return y
end
function LinearAlgebra.ldiv!(y::AbstractArray, L::Adjoint{<:Any,<:DiaSymSemiseparableCholesky}, b::AbstractArray)
    dssa_backward!(y, L.parent.Ut, L.parent.Wt, L.parent.d, b)
    return y
end
function (Base.:\)(L::DiaSymSemiseparableCholesky, b::AbstractVecOrMat)
    y = similar(b)
    dss_forward!(y, L.Ut, L.Wt, L.d, b)
    return y
end
function (Base.:\)(L::Adjoint{<:Any,<:DiaSymSemiseparableCholesky}, b::AbstractVecOrMat)
    y = similar(b)
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
		yk = @view Y[:,k]
		zk = @view Z[:,k]
		uk = @view U[:,k]
		vk = @view V[:,k]
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
    @inbounds for (u,v,w,i) in zip(eachcol(U),eachcol(V),eachcol(W),eachindex(d))
        w      .= v - P*u
        dbar[i] = sqrt(dot(u,w) + d[i])
        w      .= w/dbar[i]
        add_outer_product!(P,w) #
    end
    return W, dbar
end
#### Multiplying with Ld ####
function dss_tri_mul!(Y,U,W,ds,X)
    p     = size(U,1)
    n_rhs = size(X,2)
    Wbar  = zeros(p,n_rhs)
    @inbounds for (u,w,y,x,i) in zip(eachcol(U),eachcol(W),eachrow(Y),eachrow(X),eachindex(ds))
        add_Y_tri_diag!(y,u,Wbar,x,ds[i]) #Y[i,:] = tmpU'*Wbar + ds[i]*tmpX
        add_product!(Wbar,w,x) # Wbar  += tmpW*tmpX
    end
end
#### Adjoint of Ld ####
function dssa_tri_mul!(Y,U,W,ds,X)
    p     = size(U,1)
    n_rhs = size(X, 2)
    Ubar  = zeros(p,n_rhs)
    for (u,w,y,x,i) in zip(Iterators.reverse(eachcol(U)),
                           Iterators.reverse(eachcol(W)),
                           Iterators.reverse(eachrow(Y)),
                           Iterators.reverse(eachrow(X)),
                           Iterators.reverse(eachindex(ds)))
        add_Y_tri_diag!(y,w,Ubar,x,ds[i])    # Y[i,:] = tmpW'*Ubar + ds[i]*tmpX
        add_product!(Ubar,u,x)                    # Ubar = Ubar + tmpU*tmpX
    end
end
#### Forward substitution ####
function dss_forward!(X,U,W,ds,B)
    p     = size(U,1)
    n_rhs = size(B,2)
    Wbar  = zeros(p,n_rhs)
    for (u,w,x,b,i) in zip(eachcol(U),eachcol(W),eachrow(X),eachrow(B),eachindex(ds))
        x .= (b - Wbar'*u)/ds[i]
        add_inner_plus!(Wbar,w,x)
    end
end
#### Backward substitution ####
function dssa_backward!(X,U,W,ds,B)
    m, n = size(U)
    mx = size(B,2)
    Ubar = zeros(m,mx)
    for (u,w,b,x,i) in zip(Iterators.reverse(eachcol(U)),
                           Iterators.reverse(eachcol(W)),
                           Iterators.reverse(eachrow(B)),
                           Iterators.reverse(eachrow(X)),
                           Iterators.reverse(eachindex(ds)))
        x .= (b - Ubar'*w)/ds[i]
        add_inner_plus!(Ubar,u,x) # Ubar += tmpU .* X[i:i,:]
    end
end

#### Squared norm of columns of L = tril(UW',-1) + diag(dbar) ####
function squared_norm_cols(U,W,dbar)
    m, n = size(U)
    P = zeros(m, m)
    c = zeros(n)
    for (w,u,i) in zip(Iterators.reverse(eachcol(W)),
                       Iterators.reverse(eachcol(U)),
                       Iterators.reverse(eachindex(dbar)))
        c[i]  = dbar[i]^2 + w'*P*w
        add_outer_product!(P,u) # P += tmpU*tmpU'
    end
    return c
end
#### Implicit inverse of  L = tril(UW',-1) + diag(dbar) ####
function dss_create_yz(U, W,dbar)
    m, n = size(U)
    Y = zeros(n,m)
    Z = zeros(n,m)
    dss_forward!(Y, U, W, dbar, U')
    dssa_backward!(Z, U, W, dbar, W')
    # Probably best not to use inv
    return copy(Y'), copy((Z*inv(U*Z - Diagonal(ones(m))))')
end
