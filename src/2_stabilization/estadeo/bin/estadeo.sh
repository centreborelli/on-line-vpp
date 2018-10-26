#!/bin/bash

#this script is used to execute the program from the IPOL demo
workdir=$1
bindir=$2
tmpdir=${1}tmp
motion_method=$3
transformation=$4
smooth_strategy=$5
radius=$6
bc=$7
postprocess=$8
filetransform=transform.txt

video_in=${workdir}input_0.mp4
video_out=${workdir}output_0.mp4
graphic=${workdir}plot_params.png
file_transforms=${workdir}transform.txt
file_soft_transforms=${workdir}softtransform.txt

#create temporal directory
mkdir $tmpdir

#extract info from video
info=`avprobe -v error -show_streams  $video_in`
info="${info#*codec_type=video}"
echo

width=`echo ${info#*width=}| cut -d' ' -f 1` 
height=`echo ${info#*height=}| cut -d' ' -f 1` 
framerate=`echo ${info#*avg_frame_rate=}| cut -d' ' -f 1`
nframes=`echo ${info#*nb_frames=}| cut -d' ' -f 1`
size=${width}x${height}

if [ "$framerate" == "0/0" ] ; then
  echo "Error reading the frame rate of the video (default 30)"
  framerate="30"
fi

#show input information
echo 0.-Input parameters and video info
echo " "Smoothing radius: $radius
echo " "Dithering: $dithering
echo " "Widht: $width, Height: $height, Number of frames: $nframes, Frame rate: $framerate
echo

#convert input video to raw data rgb 3 bytes
echo 1.-Convert input video to raw data 
echo avconv -v error -i ${workdir}input_0.mp4 -f rawvideo -pix_fmt rgb24 -y $tmpdir/video.raw
time avconv -v error -i ${workdir}input_0.mp4 -f rawvideo -pix_fmt rgb24 -y $tmpdir/video.raw
echo

#execute main program
echo 2.-Estadeo video stabilization
echo $bindir/estadeo $tmpdir/video.raw $width $height $nframes -o $tmpdir/outvideo.raw -m $motion_method -t $transformation \
     -s $smooth_strategy -r $radius -b $bc -p $postprocess -w $filetransform -f $file_soft_transforms -v 
time $bindir/estadeo $tmpdir/video.raw $width $height $nframes -o $tmpdir/outvideo.raw -m $motion_method -t $transformation \
     -s $smooth_strategy -r $radius -b $bc -p $postprocess -w $filetransform -f $file_soft_transforms -v 
echo

#create output video
echo 3.-Convert output video to mp4
echo avconv -v error -f rawvideo -pix_fmt rgb24 -video_size $size -framerate $framerate -i $tmpdir/outvideo.raw  -pix_fmt yuv420p -y ${workdir}/output.mp4
time avconv -v error -f rawvideo -pix_fmt rgb24 -video_size $size -framerate $framerate -i $tmpdir/outvideo.raw  -pix_fmt yuv420p -y ${workdir}/output.mp4
echo

#create graphics
echo 4.-Generate graphics for demo
generate_graphics.sh $workdir $bindir $file_transforms $file_soft_transforms $radius

#remove temporal data
rm -R $tmpdir
