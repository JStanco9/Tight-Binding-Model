set terminal png
set output 'Bi2Se3_bands.png'
unset arrow
set style data dots
set nokey
set arrow from 0,  -1 to 0,  15.2425 nohead
set arrow from 0.468336,  -1 to 0.468336,  15.2425 nohead
set arrow from 0.643025,  -1 to 0.643025,  15.2425 nohead
set arrow from 1.11136,  -1 to 1.11136,  15.2425 nohead
set xtics (" G "  0, " L "  0.468336, " X "  0.643025, " T "  1.11136, " G"  1.28605)
set xrange [0: 1.28605]
set yrange [-1: 15.2425]
plot 'Bi2Se3_bands.dat'
set terminal xterm
plot 'Bi2Se3_bands.dat'