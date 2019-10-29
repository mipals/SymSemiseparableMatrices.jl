export SymSemiseparable

struct SymSemiseparable <: SymSemiseparableMatrix
    n::Int64
    p::Int64
    U::AbstractArray
    V::AbstractArray
end

# Constuctors
function SymSemiseparable(U::AbstractArray, V::AbstractArray)
    if size(U,1) == size(V,1) && size(U,2) == size(V,2)
        return SymSemiseparable(size(U,1),size(U,2),U,V);
    else
        error("Dimension mismatch between generators U and V")
    end
end

# Mappings
mul!(y::AbstractArray, L::SymSemiseparable, 		         x::AbstractArray) =
    ss_mul_mat!(y, L.U, L.V, x);
mul!(y::AbstractArray, L::AdjointOperator{SymSemiseparable}, x::AbstractArray) =
    ss_mul_mat!(y, L.A.U, L.A.V, x);
function inv!(y, K::SymSemiseparable, b::AbstractArray)
	L = SymSemiseparableChol(K);
	y[:,:] = L'\(L\b);
end

#####################################################################
#### Extended generator representable {p}-semiseperable matrices ####
#####################################################################

#### Matrix-matrix product ####
function ss_mul_mat!(Y::AbstractArray, U::AbstractArray, V::AbstractArray, X::AbstractArray)
    n, m = size(U);
    mx = size(X,2);
    Vbar = zeros(m,mx);
    Ubar = U'*X;
    for i = 1:n
        tmpV = V[i,:];
        tmpU = U[i,:];
        #Ubar = broadcast!(+, Ubar, tmpU .* X[i:i,:]);
        #Vbar = broadcast!(+, Vbar, tmpV .* X[i:i,:]);
        Ubar -= tmpU .* X[i:i,:];
        Vbar += tmpV .* X[i:i,:];
        Y[i,:] = Vbar'*tmpU + Ubar'*tmpV;
    end
end

#### Matrix vector product ####
function ss_create_v(U::AbstractArray, W::AbstractArray)
    n,m = size(U)
    V = zeros(n,m)
    P = zeros(m,m);
    for i = 1:n
        tmpW = W[i,:]
        tmpU = U[i,:]
        tmpV = tmpW*(tmpU'*tmpW)
        V[i,:] = tmpV + P*tmpU
        P += tmpW*tmpW';
    end
    return V
end
