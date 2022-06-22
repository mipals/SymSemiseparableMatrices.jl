"""
    alpha(p::Int)

Returns coefficients for the spline kernel of order `p`.
"""
function alpha(p::Int)
    if p < 1
        throw(DomainError("Spline order has to be strictly positive"))
    end
    a = zeros(p)
    for i = 0:p-1
        for j = 0:i
            for k = i:p-1
                a[i+1] = a[i+1] + (-1.0)^(k-j)/((k+j+1)*
                        factorial(j)*
                        factorial(p-1)*
                        factorial(p-1-k)*
                        factorial(i-j)*
                        factorial(k-i))
            end
        end
    end
    return a
end

"""
    spline_kernel(t, p)

Returns generators `U` and `V` for the spline kernel of order `p` given a list of knots `t`.
"""
function spline_kernel(t::AbstractArray, p::Int)

    fp = factorial.(p-1:-1:0)
    a  = convert.(eltype(t),alpha(p).*fp)
    Ut = (repeat(t,p,1).^Vector(p-1:-1:0))./fp
    Vt = (repeat(t,p,1).^Vector(p:2*p-1)).*a

    return Ut,Vt
end

"""
    spline_kernel_matrix(U, V)

Returns the dense spline kernel matrix with generators `U` and `V`.
"""
function spline_kernel_matrix(U, V)
    return tril(U*V') + triu(V*U',1)
end

"""
    spline_kernel_matrix(U, V, d)

Returns the dense spline kernel matrix with generators `U` and `V` plus `Diagonal(d)`.
"""
function spline_kernel_matrix(U, V, d )
    return spline_kernel_matrix(U, V) + Diagonal(d)
end

"""
    logp(Σ, T, y, σf, σn)

Generalized Log-likelihood evalation for the spline kernel.
"""
function logp(Σ, T, y, σf, σn)
    m = rank(T)
    n = length(y)
    Ki = cholesky(σf^2*Σ  + σn^2*I)
    if LinearAlgebra.issuccess(Ki) == false
        throw(DomainError("Cholesky Factorization failed"))
    end
    A = T'*(Ki\T)
    C = (Ki\T)*(A\(T'*inv(Ki)))
    lp = -0.5*y'*(Ki\y)  +
          0.5*y'*(C*y)   +
         -0.5*logdet(Ki) +
         -0.5*logdet(A)  +
         -(n-m)*0.5*log(2*π)
    return lp
end

"""
    dlogp(Σ, T, y, σf, σn)

Derivatives of the generalized Log-likelihood for the spline kernel.
"""
function dlogp(Σ, T, y, σf, σn)
    Ki = cholesky(σf^2*Σ  + σn^2*I)
    if LinearAlgebra.issuccess(Ki) == false
        throw(DomainError("Cholesky Factorization failed"))
    end
    A = T'*(Ki\T)
    # TO_DO: Solving for the identity?
    dKf = 2.0*σf*Σ
    dAf = -T'*(Ki\dKf)*(Ki\T)
    dKn = 2.0*σn*I
    dAn = -T'*(inv(Ki)*dKn)*(Ki\T)
    Pf = Ki\dKf*inv(Ki)
    Pn = Ki\(dKn*inv(Ki))
    G = T*(A\(T'*(Ki\y)))
    dlpf =  0.5*y'*(Pf*y)                +
            0.5*(-2*(y'*Pf)*G + G'*Pf*G) +
           -0.5*tr(Ki\dKf)               +
           -0.5*tr(A\dAf)
   dlpn =   0.5*y'*(Pn*y)                   +
            0.5*(-2*(y'*Pn)*G + G'*Pn*G)    +
           -0.5*tr(inv(Ki)*dKn)             +
           -0.5*tr(A\dAn)
    return [dlpf; dlpn]
end

"""
    create_T(t, m)

Creates `m`-order Taylor-polynomial basis from the knots in `t`.
"""
function create_T(t, m)
    T = ones(length(t), m)
    for ν = 2:m
        T[:, ν] = t.^(ν-1)/factorial(ν - 1)
    end
    return T
end
