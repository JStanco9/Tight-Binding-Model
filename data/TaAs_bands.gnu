set terminal png
set output 'TaAs_bands.png'
unset arrow
set style data dots
set nokey
set arrow from 0,  -1 to 0,  23.5448 nohead
set arrow from 0.525796,  -1 to 0.525796,  23.5448 nohead
set arrow from 0.674506,  -1 to 0.674506,  23.5448 nohead
set arrow from 0.823217,  -1 to 0.823217,  23.5448 nohead
set arrow from 1.2648,  -1 to 1.2648,  23.5448 nohead
set arrow from 1.55006,  -1 to 1.55006,  23.5448 nohead
set xtics (" G "  0, " S "  0.525796, " N "  0.674506, " S1 "  0.823217, " Z "  1.2648, " G "  1.55006, " X"  2.2341)
set xrange [0: 2.2341]
set yrange [-1: 23.5448]
plot 'TaAs_bands.dat'
set terminal xterm
plot 'TaAs_bands.dat'