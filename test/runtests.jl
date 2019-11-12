using SymSemiseparableMatrices
using Test
using LinearAlgebra

@testset "SymSemiseparableMatrices.jl" begin
    include("test_symsep.jl")
    include("test_symsepchol.jl")
    include("test_diasymsep.jl")
    include("test_diasymsepchol.jl")
end
