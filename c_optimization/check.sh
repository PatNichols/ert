#!/bin/bash

echo "no optimization"

make -f makefile1

time ./c_opt_0

time ./c_opt_0

time ./c_opt_1

time ./c_opt_2

time ./c_opt_2_mp

echo "some optimization"

make -f makefile2

time ./c_opt_0
time ./c_opt_1
time ./c_opt_2
time ./c_opt_2_mp

echo "max optimization"

make -f makefile3

time ./c_opt_0
time ./c_opt_1
time ./c_opt_2
time ./c_opt_2_mp

echo "profiling"

make -f makefile_profile

./c_opt_o
gprof c_opt_o >& c_opt_o.prof
./c_opt_1
gprof c_opt_1 >& c_opt_1.prof
./c_opt_2 
gprof c_opt_2 >& c_opt_2.prof