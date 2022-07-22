#!/bin/bash
make clean
make
# echo "====================================="
# ./cfd1l1.x
nohup ./cfd1l1.x &
# gnuplot script.gnu