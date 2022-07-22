#!/bin/bash
(cd ../source && make clean)
(cd ../source && make)
# echo "====================================="
(cd ../source && ./fsi3l1.x)
# nohup ./source/fsi3l1.x &
# gnuplot script.gnu