export AdjointOperator

struct AdjointOperator{T <: SymSemiseparableMatrix} <: SymSemiseparableMatrix
	A::T
	function AdjointOperator(A::T) where {T<:SymSemiseparableMatrix}
		new{T}(A)
	end
end

# Constructor
AdjointOperator(L::AdjointOperator) = AdjointOperator(L)
