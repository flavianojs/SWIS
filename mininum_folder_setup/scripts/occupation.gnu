# set terminal aqua size 450,385
# set terminal x11 size 700, 500
# set terminal canvas 
# set output "occupation.html"

set terminal png enhanced font 'Times,28' linewidth 2 size 1800, 1500
# set terminal pdf enhanced font 'Times,100' linewidth 10 size 50, 25
set output "occupation_".specifier.".png"

# set xtics ( 0.0, 8.0)
# set grid xtics

file="outputfiles/occupation_".specifier.".dat"
b=0.031
cb=1.5
b=3
a=22
a=12
# set cbrange[0.0:cb]
#et autoscale
# set xrange[-0.5:a]
# set yrange[-0.5:7]

# set xrange[-2*a/5:a]
# set yrange[-b:b]
set xrange[-a+10:a+10]
set yrange[-a:a]
set palette defined ( 1 "white", 2 "yellow", 3 "orange", 4 "red"  )
# set palette defined ( 1 "blue", 2 "green", 3 "yellow", 4 "red"  )

# plot file u 1:2:(0.46):($4)  with circles lw 0.5 lc palette fill solid 1.0 noborder notitle, "hexagonNO.dat" u 1:2 w line lw 3 notitle
# plot file u 1:2:(abs(0.46*$4)):($4)  with circles lw 0.5 lc palette fill solid 1.0 noborder notitle
plot file u 1:2:(0.25):($4)  with circles lw 0.4 lc palette fill solid 1.0 border lt 3 notitle
# plot file u 1:2:(abs(0.46*$4/GPVAL_DATA_CB_MAX)):($4)  with circles lw 0.5 lc palette fill solid 1.0 noborder notitle
