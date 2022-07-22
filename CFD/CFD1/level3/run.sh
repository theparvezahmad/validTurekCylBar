#!/bin/bash
make clean
make
# echo "====================================="
# ./cfd1l3.x
nohup ./cfd1l3.x &
# gnuplot script.gnu