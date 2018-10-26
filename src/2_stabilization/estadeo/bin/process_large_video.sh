
inc=10
second=0
for min in `seq 0 10`
do
  for sec in `seq 0 $inc 50`
  do

    start=00:$min:$sec
    ffmpeg -ss $start -i $1 -t 00:00:$inc -y tmp.mp4
    bin/cmdline_execute.sh tmp.mp4 outvideo_${min}_$sec.mp4 3 30 transforms_${min}_$sec.mat 0 "-m 0"
  done
done
