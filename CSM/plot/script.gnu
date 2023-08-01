set encoding utf8

# See https://github.com/Gnuplotting/gnuplot-palettes
# Line styles (colorbrewer Set1)
set style line 1 lc rgb '#E41A1C' pt 1 ps 1 lt 1 lw 2 # red
set style line 2 lc rgb '#377EB8' pt 6 ps 1 lt 1 lw 2 # blue
set style line 3 lc rgb '#4DAF4A' pt 2 ps 1 lt 1 lw 2 # green
set style line 4 lc rgb '#984EA3' pt 3 ps 1 lt 1 lw 2 # purple
set style line 5 lc rgb '#FF7F00' pt 4 ps 1 lt 1 lw 2 # orange
set style line 6 lc rgb '#FFFF33' pt 5 ps 1 lt 1 lw 2 # yellow
set style line 7 lc rgb '#A65628' pt 7 ps 1 lt 1 lw 2 # brown
set style line 8 lc rgb '#F781BF' pt 8 ps 1 lt 1 lw 2 # pink
# Palette
set palette maxcolors 8
set palette defined ( 0 '#E41A1C', 1 '#377EB8', 2 '#4DAF4A', 3 '#984EA3',\
4 '#FF7F00', 5 '#FFFF33', 6 '#A65628', 7 '#F781BF' )

# Standard border
set style line 11 lc rgb '#808080' lt 1 lw 3
set border 0 back ls 11
set tics out nomirror

# Standard grid
set style line 12 lc rgb '#808080' lt 0 lw 1
set grid back ls 12
unset grid

set terminal pngcairo enhanced color dashed font "Alegreya, 14" \
#rounded size 25 cm, 10 cm

# Default encoding, line styles, pallette, border and grid are set in
# /usr/local/share/gnuplot/x.y/gnuplotrc.

set xlabel "Time"
set ylabel "Deflection of Probe Point"
set grid
set key right top
#set xrange[0:6.28]
#set yrange[-1:1]
set output 'uy.png'

#plot "dynamic.dat" every ::401::2000 u 1:2 w l ls 1
#plot "dynamic.dat" every ::1601::2000 u 1:2 title "code" w l ls 1, "yDisp.txt" u 1:($2*-1) '%lf,%lf' title "benchmark" w l ls 2
#plot "dynamic.dat" every ::801::1000 u 1:2 title "code" w l ls 1, "yDisp.txt" u 1:($2*-1) '%lf,%lf' title "benchmark" w l ls 2
plot "dynamic.dat" every ::399::499 u 1:2 title "code" w l ls 1, "yDisp.txt" u 1:($2*-1) '%lf,%lf' title "benchmark" w l ls 2 dashtype 4

set output 'ux.png'
plot "dynamic.dat" every ::399::499 u 1:3 title "code" w l ls 1, "xDisp.txt" u 1:2 '%lf,%lf' title "benchmark" w l ls 2 dashtype 4
