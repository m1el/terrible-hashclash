set datafile separator ','
set term png
set logscale y
set ylabel 'time to find a collision'
set xlabel 'trail bits'
set output 'time-to-find-a-collision-singlethread.png'
plot 'data-32bits.txt' using 5 with lines title '32bits', \
     'data-40bits.txt' using 5 with lines title '40bits', \
     'data-44bits.txt' using 5 with lines title '44bits', \
     'data-48bits.txt' using 5 with lines title '48bits'
