Implementation of Local Fourier Burst Accumulation for Video Deblurring

This program is part of the IPOL publication:
    http://www.ipol.im/pub/art/2017/197/

Version 20170116

Jérémy Anger <jeremy.anger@cmla.ens-cachan.fr>, CMLA, ENS Cachan
Enric Meinhardt Llopis <enric.meinhardt@cmla.ens-cachan.fr> CMLA, ENS Cachan

This software is distributed under the terms of the BSD license (see file license.txt)

Compilation:
    run "make" to produce an executable named "fba"
    requires libpng, libtiff, libjpeg, libfftw3

Usage:
    ./fba M W P DOWNSAMPLING REGISTER ITER image1 image2 image3... REGISTERED_FORMAT OUTPUT_FORMAT

    M: half-size of the temporal window (eg: 3 to use a temporal window of 7 frames)
    W: size of patches (128 is a reasonable value, should be an even number)
    P: exponent used to compute the weight of a patch
    DOWNSAMPLING: reduction factor applied for the optical flow computation
    REGISTER: 0 to disable registration, 1 otherwise
        if the registration is disabled, REGISTERED_FORMAT is unused but should be present
    ITER: number of iterations of the algorithm
    REGISTERED_FORMAT: C-style format to indicate the path and name of registered frames (it should contain two '%d' to represent the reference frame and the registered frame)
    OUTPUT_FORMAT: C-style format to indicate the path and name of the processed frames (it should contain one '%d')

Usage example:
    (assuming a bash-like shell)
    ./fba 3 128 11 3 1 1 path/to/dataset/*.png path/to/result/registered_%03d_%03d.png path/to/result/sharp_%03d.png

    This software cannnot process video files. Please use a tool like ffmpeg to convert them (eg: ffmpeg -i input_video.avi input_image_%03d.png)

