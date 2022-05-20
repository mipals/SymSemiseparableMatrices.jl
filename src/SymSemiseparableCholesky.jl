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
                    Defining Matrix Properties / Overloading Base
==========================================================================================#
Base.Matrix(K::SymSemiseparableCholesky) = getproperty(K,:L)
Base.size(K::SymSemiseparableCholesky) = (K.n, K.n)
function Base.size(K::SymSemiseparableCholesky, d::Int)
    return (1 <= d && d <=2) ? size(K)[d] : throw(ArgumentError("Invalid dimension $d"))
end
function Base.getindex(K::SymSemiseparableCholesky, i::Int, j::Int)
	i >= j && return dot(K.Ut[:,i], K.Wt[:,j])
	return 0
end
function Base.getproperty(K::SymSemiseparableCholesky, d::Symbol)
    if d === :U
        return UpperTriangular(K.Wt'*K.Ut)
    elseif d === :L
        return LowerTriangular(K.Ut'*K.Wt)
    else
        return getfield(K, d)
    end
end
function Base.propertynames(F::SymSemiseparableCholesky, private::Bool=false)
    return (:U, :L, (private ? fieldnames(typeof(F)) : ())...)
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")},
		 K::SymSemiseparableCholesky{<:Any,<:AbstractMatrix,<:AbstractMatrix})
    summary(io, K); println(io)
    show(io, mime, K.L)
end

#==========================================================================================
                            Overloading LinearAlgebra routines
==========================================================================================#
function LinearAlgebra.det(L::SymSemiseparableCholesky)
    dd = one(eltype(L))
    @inbounds for i in 1:L.n
        dd *= dot(L.Ut[:,i],L.Wt[:,i])
    end
    return dd
end

function LinearAlgebra.logdet(L::SymSemiseparableCholesky)
    dd = zero(eltype(L))
    @inbounds for i in 1:L.n
        dd += log(dot(L.Ut[:,i],L.Wt[:,i]))
    end
    return dd
end

#==========================================================================================
                        Defining multiplication and inverse
==========================================================================================#
function LinearAlgebra.mul!(y::AbstractArray, L::SymSemiseparableCholesky, b::AbstractArray)
	ss_tri_mul!(y, L.Ut,   L.Wt,   b)
    return y
end
function LinearAlgebra.mul!(y::AbstractArray, L::Adjoint{<:Any,<:SymSemiseparableCholesky}, b::AbstractArray)
	ssa_tri_mul!(y, L.parent.Ut, L.parent.Wt, b)
    return y
end
function LinearAlgebra.ldiv!(y::AbstractArray,F::SymSemiseparableCholesky, b::AbstractArray)
	ss_forward!(y,F.Ut,F.Wt,b)
	return y
end
function LinearAlgebra.ldiv!(y::AbstractArray,F::Adjoint{<:Any,<:SymSemiseparableCholesky}, b::AbstractArray)
	ssa_backward!(y,F.parent.Ut,F.parent.Wt,b)
	return Y
end
function (Base.:\)(F::SymSemiseparableCholesky, b::AbstractVecOrMat)
    y = similar(b)
	ss_forward!(y,F.Ut,F.Wt,b)
	return y
end
function (Base.:\)(F::Adjoint{<:Any,<:SymSemiseparableCholesky}, b::AbstractVecOrMat)
    y = similar(b)
	ssa_backward!(y,F.parent.Ut,F.parent.Wt,b)
	return y
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
    @inbounds for (u,v,w) in zip(eachcol(U),eachcol(V),eachcol(W))
        w .= (v - wTw*u)/sqrt(abs(dot(u,(v - wTw*u))))
        add_outer_product!(wTw,w)
    end
    return W
end


"""
    ss_tri_mul!(Y,U,W,X)

Computes the multiplization `tril(U*W')*X = Y` in linear complexity.
"""
function ss_tri_mul!(Y,U,W,X)
    p     = size(U,1)
    n_rhs = size(X,2)
    Wbar  = zeros(p,n_rhs)
    @inbounds for (w,u,y,x) in zip(eachcol(W),eachcol(U),eachrow(Y),eachrow(X))
        add_inner_plus!(Wbar,w,x)
        add_Y_tri!(y,Wbar,u)
    end
end
"""
    ssa_tri_mul!(Y,U,W,X)

Computes the multiplization `triu(W*U')*X = Y` in linear complexity.
"""
function ssa_tri_mul!(Y,U,W,X)
    Ubar = U*X
    @inbounds for (w,u,y,x) in zip(eachcol(W),eachcol(U),eachrow(Y),eachrow(X))
        add_Y_tri!(y,Ubar,w)
        add_inner_minus!(Ubar,u,x)
    end
end
"""
    ss_forward!

Solves the linear system `tril(U*W')*X = B`.
"""
function ss_forward!(X,U,W, B)
    p     = size(U,1)
    n_rhs = size(B,2)
    Wbar  = zeros(p, n_rhs)
    @inbounds for (u,w,x,b) in zip(eachcol(U),eachcol(W),eachrow(X),eachrow(B))
        x .= (b - Wbar'*u)/dot(u,w)
        add_inner_plus!(Wbar,w,x)
    end
end
"""
    ssa_backward!(X, U, W, B)

Solves the linear system `triu(W*U')*X = Y`.
"""
function ssa_backward!(X, U, W, B)
    p     = size(U,1)
    n_rhs = size(B,2)
    Ubar  = zeros(p,n_rhs)
    @inbounds for (u,w,x,b) in zip(Iterators.reverse(eachcol(U)),
                                   Iterators.reverse(eachcol(W)),
                                   Iterators.reverse(eachrow(X)),
                                   Iterators.reverse(eachrow(B)))
        x .= (b - Ubar'w)/dot(u,w)
        add_inner_plus!(Ubar,u,x)
    end
end
