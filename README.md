# SymSemiseparableMatrices.jl

[![Build Status](https://travis-ci.com/mipals/SymSemiseparableMatrices.jl.svg?branch=master)](https://travis-ci.com/mipals/SymSemiseparableMatrices.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/mipals/SymSemiseparableMatrices.jl?svg=true)](https://ci.appveyor.com/project/mipals/SymSemiseparableMatrices-jl)
[![Codecov](https://codecov.io/gh/mipals/SymSemiseparableMatrices.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mipals/SymSemiseparableMatrices.jl)
[![Coveralls](https://coveralls.io/repos/github/mipals/SymSemiseparableMatrices.jl/badge.svg?branch=master)](https://coveralls.io/github/mipals/SymSemiseparableMatrices.jl?branch=master)

## Description
A package for efficiently computing with symmetric extended generator representable semiseparable matrices and a variant thereof. In short this means matrices of the form

<img src="https://latex.codecogs.com/svg.latex?\dpi{100}&space;K=\text{\textbf{tril}}(UV^T)&space;&plus;&space;\text{\textbf{triu}}(V^T,1),&\quad%20U,V\in\mathbb{R}^{n\times%20p}" title="K=\text{\textbf{tril}}(UV^T) + \text{\textbf{triu}}(VU^T,1)" />

as well as

<img src="https://latex.codecogs.com/svg.latex?\dpi{100}&space;K=\text{\textbf{tril}}(UV^T)&space;&plus;&space;\text{\textbf{tril}}(VU^T,1)&space;&plus;&space;\text{\textbf{diag}}(d),&\quad%20U,V\in\mathbb{R}^{n\times%20p},&space;d\in\mathbb{R}^n" title="K=\text{\textbf{tril}}(UV^T) + \text{\textbf{tril}}(VU^T,1) + \text{\textbf{diag}}(d)" />

All implemented algorithms (multiplication, Cholesky factorization, forward/backward substitution as well as various traces and determinants) run linear w.r.t. to n in time and memory using the structure of the two matrix types.

## Usage
First we need to create generators U and V that represent the symmetric matrix, ```K = tril(UV') + triu(VU',1)``` as well a test vector ```x```.
```julia
julia> U, V = spline_kernel(Vector(0.1:0.01:1), 2); # Creating input such that K is positive definite
julia> K = SymSemiseparable(U,V); # Symmetric generator representable semiseparable matrix
julia> x = ones(K.n); # Test vector
```
We can now compute products with ```K``` and ```K'```. The result are the same as ```K``` is symmetric.
```julia
julia> K*x
91×1 Array{Float64,2}:
  0.23508333333333334
  0.28261583333333334
  0.3341535          
  0.3896073333333333 
  0.44888933333333336
  0.5119124999999999 
  ⋮                  
 11.977057499999997  
 12.146079333333331  
 12.31510733333333   
 12.484138499999995  
 12.65317083333333 

julia> K'*x
91×1 Array{Float64,2}:
  0.23508333333333334
  0.28261583333333334
  0.3341535          
  0.3896073333333333 
  0.44888933333333336
  0.5119124999999999 
  ⋮                  
 11.977057499999997  
 12.146079333333331  
 12.31510733333333   
 12.484138499999995  
 12.65317083333333  

julia> K*(K\x)
91×1 Array{Float64,2}:
 1.0000000000000036
 0.9999999999999982
 0.9999999999999956
 0.9999999999999944
 0.9999999999999951
 0.999999999999995 
 ⋮                 
 0.9999999999996279
 0.9999999999996153
 0.9999999999996028
 0.9999999999995898
 0.9999999999995764
```

Furthermore from the ```SymSemiseparable``` structure we can efficiently compute the Cholesky factorization as
```julia 
julia> L = SymSemiseparableChol(K); # Computing the Cholesky factorization of K
```
Now ```L``` represents a Cholesky factorization with the form ```L = tril(UW')```. Computations with ```SymSemiseparableChol``` can be performed similar to that of ```SymSemiseparable```. The simple structure gives rise to linear (time and storage) algorithms for solving linear systems of equations.

A struct for the dealing with symmetric matrices of the form, ```K = tril(UV') + triu(VU',1) + diag(d)``` called ```DiaSymSemiseparable``` is also implemented. The usage is similar to that of ```SymSemiseparable``` and can be created as follows
```julia
julia> U, V = spline_kernel(Vector(0.1:0.01:1), 2); # Creating input such that K is positive definite
       K = DiaSymSemiseparable(U,V,rand(size(U,1)); # Symmetric generator representable semiseparable matrix + diagonal
```
The Cholesky factorization of this matrix can be computed using ```DiaSymSemiseparableChol```. Note however here that ```L``` represents a matrix of the form ```L = tril(UW',-1) + diag(c)```

## Benchmarks
### Computing Cholesky factorization of ```K = tril(UV') + triu(VU',1)```
![Scaling of the Cholesky factorization of an EGRSS matrix](https://i.imgur.com/NFqfreO.png)
### Computing Cholesky factorization of ```K = tril(UV') + triu(VU',1) + diag(d)```
![Scaling of the Cholesky factorization of an EGRQS matrix](https://i.imgur.com/IuupJSP.png)
### Solving linear systems using a Cholesky factorization with the form ```L = tril(UW')```
![Solving a system using the implicit Cholesky factorization](https://i.imgur.com/mYBNTSr.png)

## References
