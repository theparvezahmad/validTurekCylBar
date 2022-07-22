#!/bin/bash
make clean
make
echo "====================================="
# ./cfd.x
nohup ./cfd1.x &
# gnuplot script.gnu