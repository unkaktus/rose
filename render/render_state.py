#!/usr/bin/env -S python -u

import argparse

parser = argparse.ArgumentParser('render_state.py')
parser.add_argument("--total-task-number", type=int, help="Total number of tasks", default=1)
parser.add_argument("--task-id", type=int, help="Current task ID", default=0)
parser.add_argument("--state", type=str, help="State filename")
parser.add_argument("--output-dir", type=str, help="Path to output directory", default=None)
args = parser.parse_args()

import paraview.simple as pv
import numpy as np
import os

output_dir = os.path.splitext(args.state)[0]
if args.output_dir is not None:
    output_dir = args.output_dir

# Create output directory if it doesn't exist
try:
    os.mkdir(output_dir)
except:
    pass

# disable automatic camera reset on 'Show'
pv._DisableFirstRenderCameraReset()

print(f'[{args.task_id:04d}] Loading state...')
# load state
pv.LoadState(args.state)
print(f'[{args.task_id:04d}] Loaded state.')

# get animation scene
animation = pv.GetAnimationScene()


# Calculate timestamps
global_frame_times = np.linspace(animation.StartTime, animation.EndTime, animation.NumberOfFrames)
split_global_frame_times = np.array_split(global_frame_times, args.total_task_number)

frame_number_offset = 0
for a in split_global_frame_times[:args.task_id]:
    frame_number_offset += len(a)

frame_times = split_global_frame_times[args.task_id]

# Set time
for i, frame_time in enumerate(frame_times):
    global_frame_id = frame_number_offset + i
    print(f'[{args.task_id:04d}] Rendering frame #{global_frame_id:06d} ({1+i:04d} out of local batch of size {len(frame_times):04d})')
    animation.AnimationTime = frame_time
    pv.Render()
    filename = f'frame.{global_frame_id:06d}.png'
    filepath = os.path.join(output_dir, filename)
    pv.SaveScreenshot(filepath)

print(f'[{args.task_id}] Done')