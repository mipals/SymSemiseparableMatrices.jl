using SymSemiseparableMatrices
using Test
# Why do I have to import LinearAlgebra again?
using LinearAlgebra

@testset "SymSemiseparableMatrices.jl" begin
    include("test_symegrss.jl")
    include("test_symegrsscholesky.jl")
end
