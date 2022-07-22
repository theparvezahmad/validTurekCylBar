#!/bin/bash
make clean
make
# echo "====================================="
# ./cfd3.x
nohup ./cfd3l3.x &
# gnuplot script.gnu