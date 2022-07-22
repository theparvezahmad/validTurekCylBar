#!/bin/bash
make clean
make
echo "====================================="
# ./cfd.x
nohup ./cfd2.x &
# gnuplot script.gnu