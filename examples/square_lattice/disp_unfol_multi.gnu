set terminal png enhanced font 'Times,28' linewidth 2 size 900, 1400
# set output "disp_unfol_".specifier.".png"

# ymin=-0.0
# ymax=5.0
cbmax=15 #For SP
cbon=1
cbon=0

# set encoding utf8

set view map
set pm3d interpolate 0,0

set termopt enhanced
set key font "Times,20"

matriz1 = "outputfiles/disp_unfol_".specifier.".dat"
fort1 = "outputfiles/dispersion_".specifier.".dat"
fort2 = "outputfiles/dispersion_".specifier."2.dat"
fort3 = "outputfiles/const.dat"

set yrange[  ymin : ymax ]

set size 0.90, 0.35


if (cbon == 1) { set cbtics ("0" 0.0, "max" cbmax) }
if (cbon == 1) { set cbrange[-0:cbmax] }

unset xlabel
set xtics 0,.5,20 border offset 0,0.10 scale 0.0
set grid xtics lc "white" lw 0.4 lt 1

set ylabel "ω (meV)" offset 1.0,0.0
set ytics 0,1,1000  border offset 1.0, 0.0  scale 0.1 
set ytics auto  border offset 1.0, 0.0  scale 0.1 

do for [ polariz = 0:2 ] {
	if (polariz == 0) { set output "disp_unfolp1_".specifier.".png" }
	if (polariz == 1) { set output "disp_unfolp2_".specifier.".png" }
	if (polariz == 2) { set output "disp_unfolp3_".specifier.".png" }
	set multiplot

	# set colorbox user origin 1.85, 0.27 size 0.02, 0.5
	# set cbtics offset -0.8,0

	set origin 0.0, 0.66
	# set xtics ("G" 0.0, "M" 0.25, "M" 1.0, "G" 1.70710678, "Z" 2.18986540)
	load 'outputfiles/highsympoints.gnu'
	# 3 down-down, 4 down-up, 5 up-up, 6 up-down
	channel=3
	splot matriz1  u 1:2:0:(3+polariz*4) w pm3d notitle \
		, for [i=2:65] fort1 u 1:i:(0.0) w lines notitle lt 0 lc "green" lw 0.5 \

	set origin 0.0,  0.43
	# set xtics ("G" 0.0, "M" 0.25, "M" 1.0, "G" 1.70710678, "Z" 2.18986540)
	load 'outputfiles/highsympoints.gnu'
	# 3 down-down, 4 down-up, 5 up-up, 6 up-down
	channel=4
	splot matriz1  u 1:2:0:(4+polariz*4) w pm3d notitle \
		, for [i=2:65] fort1 u 1:i:(0.0) w lines notitle lt 6 lc "green" lw 1 \


	set origin 0.0,  0.2
	# set xtics ("G" 0.0, "M" 0.25, "M" 1.0, "G" 1.70710678, "Z" 2.18986540)
	load 'outputfiles/highsympoints.gnu'
	# 3 down-down, 4 down-up, 5 up-up, 6 up-down
	channel=5
	splot matriz1  u 1:2:0:5+polariz*4 w pm3d notitle \
		, for [i=2:65] fort1 u 1:i:(0.0) w lines notitle lt 1 lc "green" lw 0.5 \


	set origin 0.0, -0.03
	# set xtics ("G" 0.0, "M" 0.25, "M" 1.0, "G" 1.70710678, "Z" 2.18986540)
	load 'outputfiles/highsympoints.gnu'
	# 3 down-down, 4 down-up, 5 up-up, 6 up-down
	channel=6
	splot matriz1  u 1:2:0:6+polariz*4 w pm3d notitle \
		, for [i=2:65] fort1 u 1:i:(0.0) w lines notitle lt 1 lc "green" lw 0.5 \


	replot
	unset multiplot
} 
