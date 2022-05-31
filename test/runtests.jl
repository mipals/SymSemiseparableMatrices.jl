using Test, SafeTestsets

@safetestset "SymSemiseparableMatrix   " begin  include("test_symsep.jl")           end
@safetestset "SymSemiseparableCholesky " begin  include("test_symsepchol.jl")       end
@safetestset "DiaSymSemiseparableMatrix" begin  include("test_diasymsep.jl")        end
@safetestset "SymSemiseparableCholesky " begin  include("test_diasymsepchol.jl")    end
