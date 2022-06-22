#==========================================================================================
                                Constructors
==========================================================================================#
function DiaSymSemiseparableMatrix(U::AbstractArray, V::AbstractArray, d::AbstractArray)
    if size(U,1) == size(V,1) && size(U,2) == size(V,2) && length(d) == size(U,2)
        return DiaSymSemiseparableMatrix(size(U,2),size(U,1),U,V,d);
    else
        error("Dimension mismatch between the generators U, V and d")
    end
end
function DiaSymSemiseparableMatrix(L::SymSemiseparableMatrix, d::AbstractArray)
    return DiaSymSemiseparableMatrix(L.n, L.p, L.Ut, L.Vt, d)
end
function DiaSymSemiseparableMatrix(L::DiaSymSemiseparableCholesky)
    V, d = dss_create_vd(L.Ut, L.Wt, L.ds);
    return DiaSymSemiseparableMatrix(L.n, L.p, L.Ut, V, d)
end
#==========================================================================================
                    Defining Matrix Properties / Overloading Base
==========================================================================================#
Base.Matrix(K::DiaSymSemiseparableMatrix) = tril(K.Ut'*K.Vt) + triu(K.Vt'*K.Ut,1) + Diagonal(K.d)
Base.size(K::DiaSymSemiseparableMatrix)   = (K.n,K.n)
function getindex(K::DiaSymSemiseparableMatrix{T}, i::Int, j::Int) where T
	i > j && return dot(K.Ut[:,i],K.Vt[:,j])
	j == i && return dot(K.Vt[:,i],K.Ut[:,j]) + K.d[i]
	return dot(K.Vt[:,i],K.Ut[:,j])
end
function Base.propertynames(F::DiaSymSemiseparableMatrix, private::Bool=false)
    return (private ? fieldnames(typeof(F)) : ())
end
#==========================================================================================
                            Overloading LinearAlgebra routines
==========================================================================================#
LinearAlgebra.cholesky(K::DiaSymSemiseparableMatrix) = DiaSymSemiseparableCholesky(K)
LinearAlgebra.logdet(K::DiaSymSemiseparableMatrix)   = 2.0*logdet(cholesky(K))
LinearAlgebra.det(K::DiaSymSemiseparableMatrix)      = det(cholesky(K))^2
LinearAlgebra.logdet(L::Adjoint{<:Any,<:DiaSymSemiseparableMatrix}) = logdet(L.parent)
LinearAlgebra.det(K::Adjoint{<:Any,<:DiaSymSemiseparableMatrix}) = det(cholesky(K.parent))^2
#==========================================================================================
                        Defining multiplication and inverse
==========================================================================================#
function LinearAlgebra.mul!(y::AbstractArray, L::DiaSymSemiseparableMatrix, b::AbstractArray)
    dss_mul_mat!(y, L.Ut, L.Vt, L.d, b)
    return y
end
function LinearAlgebra.mul!(y::AbstractArray, L::Adjoint{<:Any,<:DiaSymSemiseparableMatrix}, b::AbstractArray)
    dss_mul_mat!(y, L.parent.Ut, L.parent.Vt, L.parent.d, b)
    return y
end
function LinearAlgebra.ldiv!(y::AbstractArray, K::DiaSymSemiseparableMatrix, b::AbstractArray)
	L = DiaSymSemiseparableCholesky(K)
	y .= L'\(L\b)
    return y
end
function LinearAlgebra.ldiv!(y::AbstractArray, K::Adjoint{<:Any,<:DiaSymSemiseparableMatrix}, b::AbstractArray)
	L = DiaSymSemiseparableCholesky(K.parent)
	y .= L'\(L\b)
    return y
end
function (Base.:\)(K::DiaSymSemiseparableMatrix, b::AbstractVecOrMat)
    y = similar(b)
	L = DiaSymSemiseparableCholesky(K)
	y .= L'\(L\b)
    return y
end
function (Base.:\)(K::Adjoint{<:Any,<:DiaSymSemiseparableMatrix}, b::AbstractVecOrMat)
    y = similar(b)
	L = DiaSymSemiseparableCholesky(K.parent)
	y .= L'\(L\b)
    return y
end

#===========================================================================================
                Cholesky factoriaztion of:  Higher-order quasiseparable matrices
===========================================================================================#
"""
    dss_mul_mat!(Y, U, V, d, X)

Computes the matrix-matrix product `Y = (tril(U*V') + triu(V*U',1) + diag(d))*X`.
"""
function dss_mul_mat!(Y, U, V, d, X)
    Ubar = U*X
    Vbar = zeros(eltype(X),size(Ubar))
    @inbounds for (u,v,y,x,i) in zip(eachcol(U),eachcol(V),eachrow(Y),eachrow(X),eachindex(d))
        add_inner_minus!(Ubar,u,x)
        add_inner_plus!(Vbar,v,x)
        add_Y_diag!(y,Ubar,Vbar,u,v,x,d[i]) # Y[i,:] = tmpU'*Vbar + tmpV'*Ubar + d[i]*tmpX
    end
	return Y
end

"""
    dss_create_vd(U, W, dbar)

Using `L = tril(U*W',-1) + Diagonal(dbar)`, compute `V` and `d` such that
`LL = tril(UV') + triu(V'*U,1) + Diagonal(d)`.
"""
function dss_create_vd(U, W, dbar)
    p, n = size(U)
    d = zeros(eltype(U),n)
    V = zeros(eltype(U),n,p)
    P = zeros(eltype(U),p,p)
    @inbounds for (u,w,v,i) in zip(eachcol(U),eachcol(W),eachcol(V),eachindex(dbar))
        d[i]  = dbar[i]^2 - dot(u,w)
        v    .= P*u
        add_outer_product!(P,w)
    end
    return V, d
end
