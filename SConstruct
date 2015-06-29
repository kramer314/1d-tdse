# Copyright (c) 2015 Alex Kramer <kramer.alex.kramer@gmail.com>
# See the LICENSE.txt file at the top-level directory of this distribution.

import glob

env = Environment(LINK="gfortran", F90FLAGS="-O3")
sources = glob.glob("*.f90")

env.Program("tdse.x", sources)