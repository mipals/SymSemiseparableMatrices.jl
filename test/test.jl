using SymSemiseparableMatrices
using Test
using LinearAlgebra
import SymSemiseparableMatrices: spline_kernel, spline_kernel_matrix


W = ones(3,100)
U = ones(3,100)
V = ones(3,100)
X = ones(100,2)
for (u,v,x) in zip(eachcol(U),eachcol(V),eachrow(X))
end


for (u,v) in zip(reverse(eachcol(U)),eachcol(V))
    println(dot(u,v))
end

d = rand(10)
for ds in eachrow(d)
    println(typeof(ds[1]))
end


A = rand(5,5)

for j=1:size(A,2),i = 1:size(A,1)
    println((i,j))
end
