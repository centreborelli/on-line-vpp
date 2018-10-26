#!/bin/bash

#this script is useful to execute the program from the command line
video=$1
outvideo=$2
params=$3

mkdir tmp
input_raw_video=tmp/video.raw
output_raw_video=tmp/output_video.raw

#extract info from video
info=`avprobe -v error -show_streams  $video`
info="${info#*codec_type=video}"
echo $info
echo

width=`echo ${info#*width=}| cut -d' ' -f 1` 
height=`echo ${info#*height=}| cut -d' ' -f 1` 
framerate=`echo ${info#*avg_frame_rate=}| cut -d' ' -f 1`
nframes=`echo ${info#*nb_frames=}| cut -d' ' -f 1`
size=${width}x${height}

echo Converting $video to $input_raw_video
avconv -v error -i $video -f rawvideo -pix_fmt rgb24 -y $input_raw_video
echo

echo Estadeo stabilization for video $input_raw_video to $output_raw_video
path=$(dirname "$0")
#valgrind --show-leak-kinds=all --leak-check=full --track-origins=yes $path/estadeo $input_raw_video $width $height $nframes -o $output_raw_video $params

$path/estadeo $input_raw_video $width $height $nframes -o $output_raw_video $params
echo

echo Converting $output_raw_video to $outvideo
avconv -v error -f rawvideo -pix_fmt rgb24 -video_size $size -framerate $framerate -i $output_raw_video  -pix_fmt yuv420p -y $outvideo
echo

rm -R tmp/
