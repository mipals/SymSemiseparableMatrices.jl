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

## Usage
```julia
julia> U, V = spline_kernel(Vector(0.1:0.01:1), 2); # Creating input such that K is positive definite
       K = SymSemiseparable(U,V); # Symmetric generator representable semiseparable matrix
       x = ones(K.n); # Test vector

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
 
julia> C = SymSemiseparableChol(K); # Computing the Cholesky factorization of K

julia> C*x
91×1 Array{Float64,2}:
 0.01825741858350554 
 0.022679282194091675
 0.02801063125252222 
 0.034337798127754324
 0.04166467105450973 
 0.04999152290835254 
 ⋮                   
 4.045466391180396   
 4.135793241404395   
 4.227120091628391   
 4.319446941852434   
 4.412773792076479 

julia> C'*x
91×1 Array{Float64,2}:
 12.876044456017281    
  7.289466211394393    
  4.0512284860239625   
  3.9017442701099254   
  3.8099183427776646   
  3.722847360707153    
  ⋮                    
  0.013943375672876068 
  0.009154700538161198 
  0.005366025403634955 
  0.0025773502691179007
  0.0007886751346461718

julia> C\x # Solving the linear system Cy = x using forward substitution
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

julia> C'\x # Solving the linear system C'y = x using backward substitution
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

