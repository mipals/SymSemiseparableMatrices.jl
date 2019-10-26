# SymSemiseparableMatrices.jl

[![Build Status](https://travis-ci.com/mipals/SymSemiseparableMatrices.jl.svg?branch=master)](https://travis-ci.com/mipals/SymSemiseparableMatrices.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/mipals/SymSemiseparableMatrices.jl?svg=true)](https://ci.appveyor.com/project/mipals/SymSemiseparableMatrices-jl)
[![Codecov](https://codecov.io/gh/mipals/SymSemiseparableMatrices.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mipals/SymSemiseparableMatrices.jl)
[![Coveralls](https://coveralls.io/repos/github/mipals/SymSemiseparableMatrices.jl/badge.svg?branch=master)](https://coveralls.io/github/mipals/SymSemiseparableMatrices.jl?branch=master)

A package of calculating with matrices of the form
<img src="http://www.sciweavers.org/tex2img.php?eq=K%20%3D%20%5Ctext%7B%5Ctextbf%7Btril%7D%7D%28UV%5ET%29%20%2B%20%5Ctext%7B%5Ctextbf%7Btriu%7D%7D%28VU%5ET%2C1%29&bc=White&fc=Black&im=png&fs=12&ff=arev&edit=0" align="center" border="0" alt="K = \text{\textbf{tril}}(UV^T) + \text{\textbf{triu}}(VU^T,1)" width="235" height="21" />
and
<img src="http://www.sciweavers.org/tex2img.php?eq=K%20%3D%20%5Ctext%7B%5Ctextbf%7Btril%7D%7D%28UV%5ET%29%20%2B%20%5Ctext%7B%5Ctextbf%7Btriu%7D%7D%28VU%5ET%2C1%29%20%2B%20%5Ctext%7B%5Ctextbf%7Bdiag%7D%7D%28d%29&bc=White&fc=Black&im=png&fs=12&ff=arev&edit=0" align="center" border="0" alt="K = \text{\textbf{tril}}(UV^T) + \text{\textbf{triu}}(VU^T,1) + \text{\textbf{diag}}(d)" width="321" height="22" />
