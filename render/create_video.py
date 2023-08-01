#!/usr/bin/env -S python -u
#
# Combine PNG frames into an MPEG4 video.
# Frames are selected using globbing, so missing
# frame numbers do not break the process.
#
# Optionally, in case the frames have transparent
# background, the script can set the background color.
#
# Ivan Markin, @unkaktus 2023

import os
import argparse
import subprocess
import sys
import glob

import psutil
import numpy as np
from PIL import Image


parser = argparse.ArgumentParser('create_video')
parser.add_argument("--frames-dir", type=str, help="Input directory with PNG frames", default="")
parser.add_argument("--output", type=str, help="Output file", default="output.mp4")
parser.add_argument("--frame-rate", type=int, help="Frame rate", default=10)
parser.add_argument("--background-color", type=str, help="Background color", default=None)
args = parser.parse_args()

# Use only physical cores but not more than 16
threads = int(np.min([psutil.cpu_count(logical=False), 16]))

some_frame = Image.open(glob.glob(f'{os.path.join(args.frames_dir, "*.png")}')[0])
width, height = some_frame.size

ffmpeg_cmdline = ["ffmpeg"]

if args.background_color is not None:
    ffmpeg_cmdline.extend(["-f", "lavfi",
                "-i", f"color=c={args.background_color}:s={width}x{height}"])

ffmpeg_cmdline.extend(["-vcodec", "png",
                "-framerate", f"{args.frame_rate}",
                "-pattern_type", "glob", "-i", f'{os.path.join(args.frames_dir, "*.png")}'])

if args.background_color is not None:
    ffmpeg_cmdline.extend(["-shortest", "-filter_complex", "[0:v][1:v]overlay=shortest=1,format=yuv420p[out]", "-map", "[out]"])

ffmpeg_cmdline.extend(["-pix_fmt", "yuv420p", "-crf", "17",
                "-threads", f"{threads}",
                "-preset", "fast",
                "-an",
                "-y",
                args.output])

print(' '.join(ffmpeg_cmdline))
subprocess.run(ffmpeg_cmdline, stdout=sys.stdout, stderr=sys.stderr)