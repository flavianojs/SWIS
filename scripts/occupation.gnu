# set terminal aqua size 450,385
# set terminal x11 size 700, 500
set terminal canvas 
set output "occupation.html"

# set xtics ( 0.0, 8.0)
# set grid xtics

file="occupation.dat"
a=15.5
b=0.031
b=16
a=20
# set cbrange[0.0:b]
#et autoscale
set xrange[-0.5:a]
set yrange[-0.5:7]
set xrange[-a/2:a]
set yrange[-b/2:b]
set palette defined ( 1 "blue", 2 "green", 3 "yellow", 4 "red"  )

# plot file u 1:2:(0.46):($4)  with circles lw 0.5 lc palette fill solid 1.0 noborder notitle, "hexagonNO.dat" u 1:2 w line lw 3 notitle
# plot file u 1:2:(abs(0.46*$4)):($4)  with circles lw 0.5 lc palette fill solid 1.0 noborder notitle
plot file u 1:2:(0.2):($4)  with circles lw 0.5 lc palette fill solid 1.0 noborder notitle
# plot file u 1:2:(abs(0.46*$4/GPVAL_DATA_CB_MAX)):($4)  with circles lw 0.5 lc palette fill solid 1.0 noborder notitle
