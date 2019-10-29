# SymSemiseparableMatrices.jl

[![Build Status](https://travis-ci.com/mipals/SymSemiseparableMatrices.jl.svg?branch=master)](https://travis-ci.com/mipals/SymSemiseparableMatrices.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/mipals/SymSemiseparableMatrices.jl?svg=true)](https://ci.appveyor.com/project/mipals/SymSemiseparableMatrices-jl)
[![Codecov](https://codecov.io/gh/mipals/SymSemiseparableMatrices.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mipals/SymSemiseparableMatrices.jl)
[![Coveralls](https://coveralls.io/repos/github/mipals/SymSemiseparableMatrices.jl/badge.svg?branch=master)](https://coveralls.io/github/mipals/SymSemiseparableMatrices.jl?branch=master)

## Descriptions
A package for efficiently computing with symmetric extended generator representable matrices and a variant thereof. In short this means matrices of the form

<img src="https://latex.codecogs.com/svg.latex?\dpi{100}&space;K=\text{\textbf{tril}}(UV^T)&space;&plus;&space;\text{\textbf{triu}}(VU^T,1)" title="K=\text{\textbf{tril}}(UV^T) + \text{\textbf{triu}}(VU^T,1)" />

as well as

<img src="https://latex.codecogs.com/svg.latex?\dpi{100}&space;K=\text{\textbf{tril}}(UV^T)&space;&plus;&space;\text{\textbf{tril}}(VU^T,1)&space;&plus;&space;\text{\textbf{diag}}(d)" title="K=\text{\textbf{tril}}(UV^T) + \text{\textbf{tril}}(VU^T,1) + \text{\textbf{diag}}(d)" />

All implemented algorithms (multiplication, Cholesky factorization, forward/backward substitution as well as various traces and determinants) run linear in time and memory using the structure of the two matrix types.

## Benchmarks
![Scaling of the Cholesky factorization of an EGRSS matrix](https://i.imgur.com/NFqfreO.png)
![Scaling of the Cholesky factorization of an EGRQS matrix](https://i.imgur.com/IuupJSP.png)
![Solving a system using the implicit Cholesky factorization](https://i.imgur.com/mYBNTSr.png)

## Usage
First we need to create generators U and V that represent the symmetric matrix, ```K = tril(UV') + triu(VU',1)``` as well a test vector ```x```.
```julia
julia> U, V = spline_kernel(Vector(0.1:0.01:1), 2); # Creating input such that K is positive definite
       K = SymSemiseparable(U,V); # Symmetric generator representable semiseparable matrix
       x = ones(K.n); # Test vector
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
```

Furthermore from the ```SymSemiseparable``` structure we can efficiently calculate the Cholesky factorization as
```julia 
julia> L = SymSemiseparableChol(K); # Computing the Cholesky factorization of K
```

Compuatations with ```SymSemiseparableChol``` to that of ```SymSemiseparable```, with the addition that we can solve systems
```julia 
julia> L\x # Solving the linear system Ly = x using forward substitution
91×1 Array{Float64,2}:
  54.7722557505166    
 -89.11327886790265   
  10.886919764067525  
  -2.8333003949506765 
   0.7576353467570383 
  -0.20297814193537891
   ⋮                  
   0.0                
   0.0                
   0.0                
   0.0                
   0.0    

julia> L'\x # Solving the linear system L'y = x using backward substitution
91×1 Array{Float64,2}:
   418.36207692212975   
  -425.88322666332897   
   -16.901230409582126  
    -1.1766571103485641 
    -0.0842960953287453 
    -0.00605124109635233
     ⋮                  
    30.928656802447982  
  -115.42731863148224   
   430.78061801563183   
 -1607.6951543090788    
  1267.9491923479509  

```

