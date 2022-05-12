function alpha(p::Int)
    a = zeros(p, 1)
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
    return a;
end

function spline_kernel(t::Array{Float64}, p::Int)

    #TO-DO
    # Check inputs
    # Maybe the Vector(range) can be done more eleganly.

    if all(diff(t) .> 0)
        monotonic = 1
    elseif all(diff(t) .< 0)
        monotonic = 0
    else
        throw(DomainError("t is not strictly monotonic"))
    end

    if size(t,2) != 1 && size(t,1) == 1
        t = t'
    end

    fp = factorial.(p-1:-1:0)
    a = alpha(p).*fp
    if monotonic == 1;
        U = (repeat(t,1,p).^Vector(p-1:-1:0)')./fp'
        V = (repeat(t,1,p).^Vector(p:2*p-1)').*a'
    elseif monotonic == 0
        U = (repeat(t,1,p).^Vector(p:2*p-1)')./fp'
        V = (repeat(t,1,p).^Vector(p-1:-1:0)').*a'
    end

    return U,V
end

function spline_kernel_matrix(U::Array{Float64}, V::Array{Float64})
    return tril(U*V') + triu(V*U',1)
end
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

function dlogp(Σ, T, y, σf, σn)
    m = rank(T)
    n = length(y)
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

## Auxillary functions
function create_T(t, m)
    T = ones(length(t), m)
    for ν = 2:m
        T[:, ν] = t.^(ν-1)/factorial(ν - 1)
    end
    return T
end
