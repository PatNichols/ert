#!/usr/local/bin/gnuplot -persist
# set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 600, 400 
# set output 'pm3d.11.png'
set terminal x11
set border 4095 front lt black linewidth 1.000 dashtype solid
unset parametric
set view map scale 1
set samples 50, 50
set isosamples 50, 50
unset surface
set style data lines
set xyplane relative 0
set title "nearest neighbor interpolation"
set xlabel "longitude (degress)" 
set x2range [ * : * ] noreverse writeback
set ylabel "lattitude (degrees)" 
set y2range [ * : * ] noreverse writeback
set cbrange [ * : * ] noreverse writeback
set rrange [ * : * ] noreverse writeback
set pm3d implicit at b
set palette rgbformulae 30, 31, 32
set colorbox vertical origin screen 0.9, 0.2 size screen 0.05, 0.6 front  noinvert bdefault
NO_ANIMATION = 1
splot "interp_heat.dat"
set terminal png
splot "interp_heat.dat"
