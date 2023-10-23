#!/bin/bash
#
# 1) generate grid files
#
gfortran gen_grid.f90 -o a.out && ./a.out && rm -rf a.out
#
# 2) generate xdmf file
#
gfortran gen_xdmf_flow.f90 -o a.out && ./a.out && rm -rf a.out
gfortran gen_xdmf_flow_mean.f90 -o a.out && ./a.out && rm -rf a.out
gfortran gen_xdmf_solid.f90 -o a.out && ./a.out && rm -rf a.out
