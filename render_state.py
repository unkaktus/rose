#!/usr/bin/env -S python -u

import paraview.simple as pv
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser('rose')
parser.add_argument("--total-task-number", type=int, help="Total number of tasks")
parser.add_argument("--task-id", type=int, help="Current task ID")
parser.add_argument("--state", type=str, help="State filename")
parser.add_argument("--output-dir", type=str, help="Path to output directory", default="")
args = parser.parse_args()

# disable automatic camera reset on 'Show'
pv._DisableFirstRenderCameraReset()

print(f'[{args.task_id}] Loading state...')
# load state
pv.LoadState(args.state)
print(f'[{args.task_id}] Loaded state.')

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
    print(f'[{args.task_id}] Rendering frame #{global_frame_id} ({i} out of local batch of size {len(frame_times)})')
    animation.AnimationTime = frame_time
    pv.Render()
    filename = f'frame{global_frame_id:06d}.png'
    filepath = os.path.join(args.output_dir, filename)
    pv.SaveScreenshot(filepath)

print(f'[{args.task_id}] Done')