#!/bin/bash
make clean
make
echo "====================================="
./cfdMod.x
# nohup ./cfdMod.x &
# gnuplot script.gnu