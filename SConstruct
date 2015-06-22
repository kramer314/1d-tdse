# Copyright (c) 2015 Alex Kramer <kramer.alex.kramer@gmail.com>
# See the LICENSE.txt file at the top-level directory of this distribution.

env = Environment(LINK="gfortran", LINKFLAGS="-O3")
sources = [
  "progvars.f90",
  "tridiag.f90",
  "test.f90",
]

#objs = env.Program("test.x", sources)

sources = [
  "numerics.f90",
  "progvars.f90",
  "params.f90",
  "tridiag.f90",
  "propagate.f90",
  "setup.f90",
  "tdse.f90"
]

env.Program("tdse.x", sources)