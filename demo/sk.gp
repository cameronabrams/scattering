#!/usr/bin/gnuplot
set term post eps enh col "Helvetica" 18 solid
set size square
set out "sk.eps"

set xlabel "k"
set ylabel "S(k)"

set yrange [0:4]

p "dump.sk" u 1:2 not w l
