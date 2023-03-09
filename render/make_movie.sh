#!/bin/bash
# Run inside your directory with frames named `frame.000000.png`, and this will produce a movie movie.mp4

DEFAULT_FRAME_PREFIX="frame"
FRAME_PREFIX="${1:-$DEFAULT_FRAME_PREFIX}"
OUTPUT_FILENAME="$(basename $PWD)${FRAME_PREFIX#"frame"}.mp4"

ffmpeg -vcodec png -framerate 14 -i "$FRAME_PREFIX.%06d.png" -pix_fmt yuv420p -vcodec libx264 -crf 17 -threads 16 -preset fast -y -vf "crop=trunc(iw/2)*2:trunc(ih/2)*2" $OUTPUT_FILENAME