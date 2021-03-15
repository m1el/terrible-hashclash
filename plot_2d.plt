# set datafile separator ','
set term png
# set logscale z
set xlabel 'trail bits'
set ylabel 'time to collision'
set logscale y
set xrange [0:32]
set output 'many-48.png'
plot \
     "< awk -F, 'NR>1&&$1==48{print $2,$4}' data-many.txt" with lines title '48 bits', \
     "< awk -F, 'NR>1&&$1==52{print $2,$4}' data-many.txt" with lines title '52 bits'
