import Base: adjoint, *, \
import LinearAlgebra: logdet

###### ' ######
adjoint(L::T) where {T <: SymSemiseparableMatrix} = AdjointOperator(L)

## Fix (*) and (\) overloadings. Not nice to have if/else statements
###### * ######
function (*)(L::SymSemiseparableMatrix, x::AbstractArray)
	if typeof(L) <: AdjointOperator
		y = zeros(L.A.n, size(x,2));
		mul!(y, L, x)
		return y
	else
		y = zeros(L.n, size(x,2));
		mul!(y, L, x)
		return y
	end
end

###### \ ######
function (\)(L::SymSemiseparableMatrix, x::AbstractArray)
	if typeof(L) <: AdjointOperator
		y = zeros(L.A.n, size(x,2));
		inv!(y, L, x)
		return y
	else
		y = zeros(L.n, size(x,2));
		inv!(y, L, x)
		return y
	end
end
#### Log-determinant ####
function logdet(L::SymSemiseparableMatrix)
	return newlogdet(L)
end
