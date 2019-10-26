# SymSemiseparableMatrices.jl

[![Build Status](https://travis-ci.com/mipals/SymSemiseparableMatrices.jl.svg?branch=master)](https://travis-ci.com/mipals/SymSemiseparableMatrices.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/mipals/SymSemiseparableMatrices.jl?svg=true)](https://ci.appveyor.com/project/mipals/SymSemiseparableMatrices-jl)
[![Codecov](https://codecov.io/gh/mipals/SymSemiseparableMatrices.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mipals/SymSemiseparableMatrices.jl)
[![Coveralls](https://coveralls.io/repos/github/mipals/SymSemiseparableMatrices.jl/badge.svg?branch=master)](https://coveralls.io/github/mipals/SymSemiseparableMatrices.jl?branch=master)

A package for efficiently calculating with matrices of the form
<a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{100}&space;K=\text{\textbf{tril}}(UV^T)&space;&plus;&space;\text{\textbf{triu}}(VU^T,1)" target="_blank"><img src="https://latex.codecogs.com/png.latex?\dpi{100}&space;K=\text{\textbf{tril}}(UV^T)&space;&plus;&space;\text{\textbf{triu}}(VU^T,1)" title="K=\text{\textbf{tril}}(UV^T) + \text{\textbf{triu}}(VU^T,1)" /></a>
and 
<a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{100}&space;K=\text{\textbf{tril}}(UV^T)&space;&plus;&space;\text{\textbf{tril}}(VU^T,1)&space;&plus;&space;\text{\textbf{diag}}(d)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{100}&space;K=\text{\textbf{tril}}(UV^T)&space;&plus;&space;\text{\textbf{tril}}(VU^T,1)&space;&plus;&space;\text{\textbf{diag}}(d)" title="K=\text{\textbf{tril}}(UV^T) + \text{\textbf{tril}}(VU^T,1) + \text{\textbf{diag}}(d)" /></a>
All algorithsm run linear in time and memory.
