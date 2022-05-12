# Overloading adjoint
adjoint(L::T) where {T <: SymSemiseparableMatrix} = AdjointOperator(L)

# TODO: Fix (*) and (\) overloadings. Not nice to have if/else statements

# Overloading multiplication
function (*)(L::SymSemiseparableMatrix, x::AbstractArray)
	if typeof(L) <: AdjointOperator
		y = zeros(L.A.n, size(x,2));
		mul!(y, L, x)
	else
		y = zeros(L.n, size(x,2));
		mul!(y, L, x)
	end
    return y
end

# Overloading inverse
function (\)(L::SymSemiseparableMatrix, x::AbstractArray)
	if typeof(L) <: AdjointOperator
		y = zeros(L.A.n, size(x,2));
		inv!(y, L, x)
	else
		y = zeros(L.n, size(x,2));
		inv!(y, L, x)
	end
    return y
end
# Overloading log-determinant
function logdet(L::SymSemiseparableMatrix)
	return newlogdet(L)
end
