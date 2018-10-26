#!/bin/bash

dir=$1
bindir=$2
transforms=$3
softtransforms=$4
radius=$5

#------------ GENERTATE DATA FOR GRAPHICS -------------
#generate velocity data
generate_graphics v $transforms $dir/velocity.txt
generate_graphics v $softtransforms $dir/softvelocity.txt

#generate rotation and zoom evolution data
generate_graphics r $transforms $dir/rotation.txt
generate_graphics z $transforms $dir/zoom.txt

#generate smooth trajectories from points
generate_graphics s $dir/velocity.txt $dir/softpoints.txt $radius
generate_graphics m $dir/zoom.txt $dir/softzoomscalar.txt $radius
generate_graphics m $dir/rotation.txt $dir/softrotationscalar.txt $radius

#generate linear scale-space from smooth points
generate_graphics l $dir/softpoints.txt $dir/velocity.txt $dir/linear_scale_space.txt
generate_graphics e $dir/softzoomscalar.txt $dir/zoom.txt $dir/linear_scale_space_zoom.txt
generate_graphics e $dir/softrotationscalar.txt $dir/rotation.txt $dir/linear_scale_space_rotation.txt

#generate dft data from linear scale space
cut -d " " -f 1 $dir/linear_scale_space.txt > .tmp.txt
generate_graphics f .tmp.txt $dir/linear_scale_spacex.dft
cut -d " " -f 2 $dir/linear_scale_space.txt > .tmp.txt
generate_graphics f .tmp.txt $dir/linear_scale_spacey.dft

generate_graphics f $dir/linear_scale_space_rotation.txt $dir/linear_scale_space_rotation.dft
generate_graphics f $dir/linear_scale_space_zoom.txt $dir/linear_scale_space_zoom.dft

rm .tmp.txt

#---------------- GENERATE GRAPHICS ------------------
gnuplot -e "file='${dir}velocity.png'; file0='${dir}velocity.txt'; \
file1='${dir}softvelocity.txt'" ${bindir}plgs/velocity.plg

#Generate linear scale space graphics
gnuplot -e "file='${dir}linear_scale_space.png'; file0='${dir}velocity.txt'; \
file1='${dir}linear_scale_space.txt'; file2='${dir}softpoints.txt'" ${bindir}plgs/scale_space2.plg

#Generate zoom linear scale space graphics
gnuplot -e "file='${dir}zoom_linear_scale_space.png'; file0='${dir}zoom.txt'; \
file1='${dir}linear_scale_space_zoom.txt'; file2='${dir}softzoomscalar.txt'; \
titulo='z'; constant=1" ${bindir}plgs/scale_space.plg

#Generate rotation linear scale space graphics
gnuplot -e "file='${dir}rotation_linear_scale_space.png'; file0='${dir}rotation.txt'; \
file1='${dir}linear_scale_space_rotation.txt'; file2='${dir}softrotationscalar.txt'; \
titulo='angle'; constant=0" ${bindir}plgs/scale_space.plg

#Generate linear dft graphics
gnuplot -e "file='${dir}linear_dft.png'; file10='${dir}linear_scale_spacex.dft'; \
file11='${dir}linear_scale_spacey.dft'" ${bindir}plgs/dft2.plg

#Generate zoom linear dft graphics
gnuplot -e "file='${dir}zoom_linear_dft.png'; \
file1='${dir}linear_scale_space_zoom.dft'; titulo='z'" ${bindir}plgs/dft.plg

#Generate rotation linear dft graphics
gnuplot -e "file='${dir}rotation_linear_dft.png'; \
file1='${dir}linear_scale_space_rotation.dft'; titulo='angle" ${bindir}plgs/dft.plg

#rm $dir/*.txt*
rm $dir/*.dft