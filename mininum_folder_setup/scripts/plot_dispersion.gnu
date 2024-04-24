set terminal png enhanced font 'Times,28' linewidth 2 size 900, 500
# set output "disp_unfol_".specifier.".png"

# ymin=-0.0
# ymax=5.0
cbmax=15 #For SP
cbmax=5 #For SP
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

# set size 1.02, 0.80
# set size 1.02, 0.40

# set size 0.34, 0.35
# set size 1.04, 0.35
# set size 0.92, 0.45
set size 1.15, 1.25


if (cbon == 1) { set cbtics ("0" 0.0, "max" cbmax) }
if (cbon == 1) { set cbrange[-0:cbmax] }

unset xlabel
set xtics 0,.5,20 border offset 0,0.10 scale 0.0
set grid xtics lc "white" lw 0.4 lt 1

set ylabel "ω (meV)" offset 1.0,0.0
set ytics 0,125,1000 border offset 1.0, 0.0  scale 0.1 
set ytics 0,200,1000  border offset 1.0, 0.0  scale 0.1 
set ytics 0,10,1000  border offset 1.0, 0.0  scale 0.1 
set ytics 0,1,1000  border offset 1.0, 0.0  scale 0.1 
set ytics auto  border offset 1.0, 0.0  scale 0.1 
# set ytics 0,0.5,1000 border offset 1.0, 0.0  scale 0.1 
# set grid ytics lc "white" lw 0.4 lt 1

do for [ polariz = 1:1 ] {
	if (polariz == 0) { set output "disp_unfolpX_".specifier.".png" }
	if (polariz == 1) { set output "disp_unfolpY_".specifier.".png" }
	if (polariz == 2) { set output "disp_unfolpZ_".specifier.".png" }
	set multiplot

	aa=4+polariz*4
	bb=6+polariz*4

	set colorbox user origin 1.85, 0.27 size 0.02, 0.5
	set cbtics offset -0.8,0

   # set xrange[0:1]
	set yrange[  ymin : 9 ]
	set yrange[  ymin : ymax ]
	# set origin 0.03, 0.60
	set origin -0.01, -0.15
	set xtics ("L" 0.0, "Γ" 0.86602540, "X" 1.86602540, "W" 2.36602540, "Γ" 3.48405939)
	load 'outputfiles/highsympoints.gnu'
	# 3 down-down, 4 down-up, 5 up-up, 6 up-down
	n_a1 = 20
	natom_small_cell = 3
	naucell=natom_small_cell*n_a1
	mode_a=1       #1  2  3
	mode_b=n_a1+1+4   #21 22
	mode_c=2*n_a1+3   #40 41 42
	# mode_b=n_a1+1   #21 22
	# mode_c=2*n_a1   #40 41 42
	# mode_b=3*n_a1+1   #21 22
	# mode_c=4*n_a1   #40 41 42
	mode_d=naucell #60
	splot matriz1  u 1:2:0:($7+$8+$9+$10) w pm3d notitle \
		, for [i=2:300] fort1 u 1:i:(0.0) w lines notitle lt 6 lc "green" lw 0.1  \
		# , for [i=2:300] fort1 u 1:i:(0.0) w lines notitle lt 6 lc "green" lw 0.1  \
		, for [i=(naucell+2-2)-mode_c:(naucell+2)-mode_c] fort1 u 1:i:(0.0) w lines notitle lt 6 lc "blue" lw 0.5  \
		, for [i=(naucell+2-1)-mode_b:(naucell+2)-mode_b] fort1 u 1:i:(0.0) w lines notitle lt 6 lc "blue" lw 0.5  \
		# , for [i=(naucell+2-2)-mode_a:(naucell+2)-mode_a] fort1 u 1:i:(0.0) w lines notitle lt 6 lc "blue" lw 0.5  \
		# , for [i=(naucell+2-0)-mode_d:(naucell+2)-mode_d] fort1 u 1:i:(0.0) w lines notitle lt 6 lc "blue" lw 0.5  \

	replot
	unset multiplot
} 
