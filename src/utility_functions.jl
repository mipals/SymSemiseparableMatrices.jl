"""
add_Y_diag!(Y,Ubar,Vbar,u,v,x,d)

Inplace computation of `Y[i,:] = Vbar'*u + Ubar'*v + d[i]*x`.
"""
function add_Y_diag!(Y,Ubar,Vbar,u,v,x,d)
    @inbounds for k = 1:size(Vbar,2)
        Y[k] = d*x[k]
        @inbounds for i = 1:size(Vbar,1)
            Y[k] += Vbar[i,k]*u[i] + Ubar[i,k]*v[i]
        end
    end
end
"""
add_Y_tril_diag!(Y,w,Ubar,x,d)

Inplace computation of `Y[i,:] = w'*Ubar + ds[i]*x`.
"""
function add_Y_tri_diag!(Y,w,Ubar,x,d)
    @inbounds for k = 1:size(Ubar,2)
        Y[k] = d*x[k]
        @inbounds for i = 1:size(Ubar,1)
            Y[k] += Ubar[i,k]*w[i]
        end
    end
end
"""
add_Y_tri!(Y,Wbar,u)

Inplace computation of `Y[i,:] = Wbar'*u`.
"""
function add_Y_tri!(Y,Wbar,u)
    for k = 1:size(Wbar,2)
        Y[k] = 0.0
        for i = 1:size(Wbar,1)
            Y[k] += Wbar[i,k]*u[i]
        end
    end
end

"""
add_Y!(Y,Ubar,Vbar,u,v)

Inplace computation of `Y[i,:] = Vbar'*u + Ubar'*v`.
"""
function add_Y!(Y,Ubar,Vbar,u,v)
    @inbounds for k = 1:size(Vbar,2)
        Y[k] = 0.0
        @inbounds for i = 1:size(Vbar,1)
            Y[k] += Vbar[i,k]*u[i] + Ubar[i,k]*v[i]
        end
    end
end

"""
add_inner_minus!(Ubar,u,x)

Inplace computation of `Ubar -= u .* x`.
"""
function add_inner_minus!(Ubar,u,x)
    for j=1:size(Ubar,2),i = 1:size(Ubar,1)
        Ubar[i,j] = Ubar[i,j] - u[i]*x[j]
    end
end
"""
add_inner_plus!(Vbar,v,x)

Inplace computation of `Vbar += v .* x`.
"""
function add_inner_plus!(Vbar,v,x)
    for j=1:size(Vbar,2),i = 1:size(Vbar,1)
        Vbar[i,j] = Vbar[i,j] + v[i]*x[j]
    end
end
"""
    add_outer_product!(wTw,w)

Inplace computation of `wTw = wTw + w*w'`
"""
function add_outer_product!(wTw,w)
    for j=1:size(wTw,2),i = 1:size(wTw,1)
        wTw[i,j] = wTw[i,j] + w[i]*w[j]
    end
end
"""
add_product!(Wbar,w,x)

Inplace computation of `Wbar  += w*x'`
"""
function add_product!(Wbar,w,x)
    for j=1:size(Wbar,2),i = 1:size(Wbar,1)
        Wbar[i,j] = Wbar[i,j] + w[i]*x[j]
    end
end
