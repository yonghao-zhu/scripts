#!/bin/sh
#edit by lipai@USTC
#plot energies of structure in optimization
#awk '/E0/{if ( i<=5 ) i++;else print $0 }' OSZICAR >temp.e
gnuplot <<EOF
set term dumb
set title 'E_tot-Time'
set xlabel 'Time (fs)'
set ylabel 'E_tot(eV)'
plot 'MDSTEPS' u 2:4 w l
EOF

