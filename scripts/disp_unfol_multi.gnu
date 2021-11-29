set terminal png enhanced font 'Times,26' linewidth 2 size 1000, 1000

cbmax=0.20
cbon=0

# set encoding utf8

set view map
set pm3d interpolate 0,0

set termopt enhanced
set key font "Times,18"

matriz1 = "outputfiles/disp_unfol_".specifier.".dat"
fort1 = "outputfiles/dispersion_".specifier.".dat"
fort2 = "outputfiles/dispersion_".specifier."2.dat"
fort3 = "outputfiles/const.dat"

set yrange[  ymin : ymax ]

# set size 1.1, 0.7
set size 1.1, 1.1
set size 1.07, 0.5



if (cbon == 1) { set cbtics ("0" 0.0, "max" cbmax) }
if (cbon == 1) { set cbrange[-0:cbmax] }

unset xlabel
set xtics 0,.5,20 border offset 0,0.10 scale 0.0
set grid xtics lc "white" lw 0.4 lt 1

set ylabel "ω (meV)" offset 1.5,0.0
set ytics 0,125,1000 border offset 1.0, 0.0  scale 0.1 
set ytics 0,200,1000  border offset 1.0, 0.0  scale 0.1 
set ytics 0,1,1000  border offset 1.0, 0.0  scale 0.1 
set ytics 0,10,1000  border offset 1.5, 0.0  scale 0.1 
# set ytics 0,1,1000  border offset 1.0, 0.0  scale 0.1 
# set ytics 0,0.5,1000 border offset 1.0, 0.0  scale 0.1 
# set ytics 0,0.1,1000 border offset 1.0, 0.0  scale 0.1 
# set grid ytics lc "white" lw 0.4 lt 1

do for [ polariz = 2:2 ] {
	if (polariz == 0) { set output "disp_unfolp1_".specifier.".png" }
	if (polariz == 1) { set output "disp_unfolp2_".specifier.".png" }
	if (polariz == 2) { set output "disp_unfolp3_".specifier.".png" }
	set multiplot

	aa=4+polariz*4
	bb=6+polariz*4

	set colorbox user origin 0.90, 0.27 size 0.02, 0.5
	# set colorbox size 0.02, 0.5
	set cbtics offset -0.8,0
	
	set origin -0.05, 0.56
	# set origin -0.07, 0.16
	# set origin -0.07, 0.38
	load 'outputfiles/highsympoints.gnu'
	# set xtics ("M"  0.0, "Γ"  0.08455945, "K"  0.18219936, "M"  0.23101795, "Γ"  0.31557740, "Z"  0.38655588)
	# 3 down-down, 4 down-up, 5 up-up, 6 up-down
	channel=3
	splot matriz1  u 1:2:0:channel+polariz*4 w pm3d notitle \
		# , for [i=2:1000] fort1 u 1:i:(0.0) w lines notitle lt 1 lc "green" lw 0.5 \

	# unset colorbox

	set origin -0.05, 0.25
	# set xtics ("G" 0.0, "M" 0.25, "M" 1.0, "G" 1.70710678, "Z" 2.18986540)
	load 'outputfiles/highsympoints.gnu'
	# 3 down-down, 4 down-up, 5 up-up, 6 up-down
	channel=4
	splot matriz1  u 1:2:0:channel+polariz*4 w pm3d notitle \
		# , for [i=2:1000] fort1 u 1:i:(0.0) w lines notitle lt 1 lc "green" lw 0.5 \
		# , fort3 u 1:2:0 w lines notitle lc "white" \


	set origin -0.05, -0.06
	# set xtics ("G" 0.0, "M" 0.25, "M" 1.0, "G" 1.70710678, "Z" 2.18986540)
	load 'outputfiles/highsympoints.gnu'
	set xtics ("M"  0.0, "G"  0.07216878, "M"  0.14193194)
	# 3 down-down, 4 down-up, 5 up-up, 6 up-down
	channel=6
	splot matriz1  u 1:2:0:channel+polariz*4 w pm3d notitle \
		# , for [i=2:1000] fort1 u 1:i:(0.0) w lines notitle lt 1 lc "green" lw 0.5 \
		# , fort3 u 1:2:0 w lines notitle lc "white" \

	replot
	unset multiplot
} 
