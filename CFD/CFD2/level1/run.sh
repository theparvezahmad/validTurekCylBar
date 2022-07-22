#!/bin/bash
make clean
make
# echo "====================================="
# ./cfd.x
nohup ./cfd2l1.x &
# gnuplot script.gnu