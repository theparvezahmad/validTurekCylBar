#!/bin/bash
make clean
make
echo "====================================="
# ./cfd.x
nohup ./cfd2l3.x &
# gnuplot script.gnu