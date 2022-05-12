using SymSemiseparableMatrices
using Test
using LinearAlgebra

import SymSemiseparableMatrices: spline_kernel, spline_kernel_matrix

@testset "SymSemiseparableMatrices.jl" begin
    include("test_symsep.jl")
    include("test_symsepchol.jl")
    include("test_diasymsep.jl")
    include("test_diasymsepchol.jl")
end
