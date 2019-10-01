# Video Restoration Pipeline

This an *on-line* video processing pipeline with the objective of enhancing the quality of a video, built around [vpp](https://github.com/kidanger/vpp). By on-line we mean that it works a frame-by-frame. It consists of a series of processes connected by unix pipes. These are some of the available tools:
- video stabilization
- optical flow
- noise estimation
- band noise removal
- fixed pattern noise removal
- AWGN denoising
- tone mapping


The following pipeline is an example:
```
mkfifo noisy
mkfifo oflow
sigma=20
bin/readvid "/path/to/vid/*png" - \
	| bin/vlambda - "randn $sigma * x +" \
	| bin/vp dup - noisy - \
	| bin/tvl1flow - - \
	| bin/vp dup - oflow - \
	| bin/vlambda - "x(0,0)[0] x(-1,0)[0] - x(0,0)[1] x(0,-1)[1] - + 0.5 > 255 *" -o - \
	| bin/kalman -i noisy -o oflow -k - -s $sigma -d - \
	| bin/writevid - /tmp/d%03d.png
rm noisy oflow
```

It creates the following pipeline:
```
readvid 🠆 add noise 🠆 dup 🠆 tvl1flow 🠆 dup 🠆 occlusions 🠆 kalman 🠆 writevid
                       │                │                   🡑  🡑
                       │                └────(oflow)────────┘  │
                       └─────(noisy)───────────────────────────┘
```






