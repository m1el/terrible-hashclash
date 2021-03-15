set datafile separator ','
set term png
set xlabel 'trail bits'
set ylabel 'time to collision'
set logscale y
set xrange [48:64]
set output 'best-hpc.png'
plot 'best-hashes-per-coll.txt' using 1:3 with lines title 'best hashes per collision'
