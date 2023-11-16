#!/bin/sh
#edit by lipai@USTC
#plot energies of structure in optimization
#awk '/E0/{if ( i<=5 ) i++;else print $0 }' OSZICAR >temp.e
gnuplot <<EOF
set term dumb
set title 'Temp-Time'
set xlabel 'Time (fs)'
set ylabel 'TEMP(K)'
plot 'MDSTEPS' u 2:8 w l
EOF

