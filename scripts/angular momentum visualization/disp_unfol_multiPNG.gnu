set terminal png enhanced font 'Times,28' linewidth 2 size 600,600
set output "disp_unfol_".specifier.".png"

# ymin=-0.0
# ymax=5.0
cbmax=500 #For AS
cbmax=1000 #For SsK
cbmax=22 #For SsK
cbon=0

set encoding utf8

set view map
set pm3d interpolate 0,0

set termopt enhanced
set key font "Times,20"

matriz1 = "outputfiles/disp_unfol_".specifier.".dat"
fort1 = "outputfiles/dispersion_".specifier.".dat"
matriz2 = "outputfiles/disp_unfol_".specifier.".dat"

set cbtics offset -0.8,0

set autoscale
# set xrange[0:1.] 
set yrange[  ymin : ymax ]
set xtics 0,.5,20
set xtics border offset 0,0.5
set xtics ("Γ" 0.0, "X" 0.5, "Γ" 0.99, "M" 1.5, "Γ" 1.9875)
if (cbon == 1) { set cbtics ("0" 0.0, "max" cbmax) }
set grid xtics lc "white" lw 0.4 lt 1
set xtics scale 0.0
if (cbon == 1) { set cbrange[-0:cbmax] }

# set palette defined ( 1 "white",2 "yellow",3 "green", 4 "blue" )

set multiplot

# set size 0.57, 1.3

set ytics border scale 0.5

# unset ylabel
# unset ytics
# set origin 0.43, -0.15
# set xtics ("X" 0.0, "M" 0.5, "X" 0.99375000)
# splot matriz2  u 1:2:0:($3+$4+$5+$6) w pm3d notitle \

# unset colorbox
set ylabel "ω/J" offset 0.,0.0
set ytics scale 0.5 offset 0.3, 0.0
set ytics 0,2,20
# set origin -0.01, -0.15
set xtics ("Y" 0.0, "M" 0.5, "Y" 0.99375000)

# splot matriz1  u 1:2:0:($7) w pm3d notitle \

splot matriz1  u 1:2:0:($3+$4+$5+$6) w pm3d notitle \
   , for [i=2:65] fort1 u 1:i:(0.0) w lines notitle lt 2 lw 0.5
