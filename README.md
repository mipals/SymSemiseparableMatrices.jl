# SymSemiseparableMatrices.jl

[![CI](https://github.com/mipals/SymSemiseparableMatrices.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/mipals/SymSemiseparableMatrices.jl/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/mipals/SymSemiseparableMatrices.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mipals/SymSemiseparableMatrices.jl)

## Description
A package for efficiently computing with symmetric extended generator representable semiseparable matrices (EGRSS Matrices) and a varient with added diagonal terms. In short this means matrices of the form
$$
K = \text{\textbf{tril}}(UV^\top) + \text{\textbf{triu}}\left(VU^\top,1\right), \quad U,V\in\mathbb{R}^{p\times n}
$$

$$
K = \text{\textbf{tril}}(UV^\top) + \text{\textbf{triu}}\left(VU^\top,1\right) + \text{\textbf{diag}}(d), \quad U,V\in\mathbb{R}^{p\times n},\ d\in\mathbb{R}^n
$$

All implemented algorithms (multiplication, Cholesky factorization, forward/backward substitution as well as various traces and determinants) scales with $O(p^kn)$. Since $p \ll n$ this result in very scalable computations.

A more in-depth descriptions of the algorithms can be found in [1] or [here](https://github.com/mipals/SmoothingSplines.jl/blob/master/mt_mikkel_paltorp.pdf).

## Usage
Adding the package can be done through
```
(@v1.7) pkg> add SymSemiseparableMatrices
```
First we need to create generators U and V that represent the symmetric matrix, ```K = tril(UV') + triu(VU',1)``` as well a test vector ```x```.
```julia
julia> using SymSemiseparableMatrices
julia> import SymSemiseparableMatrices: spline_kernel
julia> U, V = spline_kernel(Vector(0.1:0.01:1)', 2); # Creating input such that K is PD
julia> K = SymSemiseparableMatrix(U,V); # Symmetric generator representable semiseparable matrix
julia> x = ones(size(K)); # Test vector
```
We can now compute products with ```K``` and ```K'```. The result are the same, since ```K``` is symmetric.
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

Furthermore from the ```SymSemiseparableMatrix``` structure we can efficiently compute the Cholesky factorization as
```julia 
julia> L = cholesky(K); # Computing the Cholesky factorization of K
julia> K*(L'\(L\x))
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
Now ```L``` represents a Cholesky factorization with of form ```L = tril(UW')```, requiring only $O(pn)$ storage. 

A struct for the dealing with symmetric matrices of the form, ```K = tril(UV') + triu(VU',1) + diag(d)``` called ```DiaSymSemiseparableMatrix``` is also implemented. The usage is similar to that of ```DiaSymSemiseparableMatrix``` and can be created as follows
```julia
julia> U, V = spline_kernel(Vector(0.1:0.01:1)', 2); # Creating input such that K is PD
julia> K = DiaSymSemiseparableMatrix(U,V,rand(size(U,2)); # Symmetric EGRSS matrix + diagonal
```
The Cholesky factorization of this matrix can be computed using ```cholesky```. Note however here that ```L``` represents a matrix of the form ```L = tril(UW',-1) + diag(c)```

## Benchmarks
### Computing Cholesky factorization of ```K = tril(UV') + triu(VU',1)```
![Scaling of the Cholesky factorization of an SymSemiseparableMatrix matrix](https://i.imgur.com/NFqfreO.png)
### Computing Cholesky factorization of ```K = tril(UV') + triu(VU',1) + diag(d)```
![Scaling of the Cholesky factorization of an DiaSymSemiseparableMatrix matrix](https://i.imgur.com/IuupJSP.png)
### Solving linear systems using a Cholesky factorization with the form ```L = tril(UW')```
![Solving a system using the implicit Cholesky factorization](https://i.imgur.com/mYBNTSr.png)

## References
[1] M. S. Andersen and T. Chen, “Smoothing Splines and Rank Structured Matrices: Revisiting the Spline Kernel,” SIAM Journal on Matrix Analysis and Applications, 2020.

[2] J. Keiner. "Fast Polynomial Transforms." Logos Verlag Berlin, 2011.
