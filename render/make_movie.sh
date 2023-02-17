#!/bin/bash
# Run inside your directory with frames named `frame.000000.png`, and this will produce a movie movie.mp4

ffmpeg -vcodec png -framerate 14 -i frame.%06d.png -pix_fmt yuv420p -vcodec libx264 -crf 17 -threads 64 -preset fast -y -vf "crop=trunc(iw/2)*2:trunc(ih/2)*2" movie.mp4