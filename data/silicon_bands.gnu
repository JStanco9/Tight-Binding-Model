set terminal png
set output 'silicon_bands.png'
unset arrow
set style data dots
set nokey
set arrow from 0,  -6.6983 to 0,  17.6539 nohead
set arrow from 0.53347,  -6.6983 to 0.53347,  17.6539 nohead
set arrow from 1.14947,  -6.6983 to 1.14947,  17.6539 nohead
set arrow from 1.36726,  -6.6983 to 1.36726,  17.6539 nohead
set xtics (" L "  0, " G "  0.53347, " X "  1.14947, " K "  1.36726, " G"  2.02062)
set xrange [0: 2.02062]
set yrange [-6.6983: 17.6539]
plot 'silicon_bands.dat'
set terminal xterm
plot 'silicon_bands.dat'