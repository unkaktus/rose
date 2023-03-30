#!/usr/bin/env -S python -u

import argparse

parser = argparse.ArgumentParser('overploy.py')
parser.add_argument("--total-task-number", type=int, help="Total number of tasks", default=1)
parser.add_argument("--task-id", type=int, help="Current task ID", default=0)
parser.add_argument("--state", type=str, help="State filename")
parser.add_argument("--output-dir", type=str, help="Path to output directory", default=None)

parser.add_argument("--frame-spacing", type=float, help="Frame spacing in total masses", default=10)

parser.add_argument("--m1", type=float, help="Mass of the primary component", default=None)
parser.add_argument("--m2", type=float, help="Mass of the secondary component", default=None)

# TODO: Can be taken on from the waveform
parser.add_argument("--merger-time", type=float, help="Merger time", default=None)

parser.add_argument("--astro-units-flux", type=bool, help="Use astrophysical units for flux (L_sun/sr)", default=False)

args = parser.parse_args()

import os
import math
import numpy as np
import paraview.simple as pv

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FuncFormatter
from astropy import units as u, constants as const


# better pictures and legends
plt.rc('figure', dpi=250)
plt.rc('text', usetex=True)

class GU():
    def __init__(self, v):
        self.v = v*u.Msun
        self.M = self._to(u.g)
        self.L = self._to(u.cm)
        self.T = self._to(u.s)
        self.Density = self.MLT(1, -3, 0)
        self.Energy = self.MLT(1, 2, -2)
        self.MagneticFluxDensity = self.MLT(1/2, -1/2, -1)
        self.Luminosity = self.Energy/self.T
    def _to(self, unit):
        if unit.is_equivalent(u.Msun):
            return self.v.to(unit)
        if unit.is_equivalent(u.cm):
            return (self.v*const.G/const.c**2).to(unit)
        if unit.is_equivalent(u.s):
            return (self.v*const.G/const.c**3).to(unit)
    def MLT(self, m, l, t):
        return self.M**m * self.L**l * self.T**t
    

def time_text(time, mass, merger) -> str:
    t = GU(mass).T.to('ms').value * (time - merger)
    if t < 0:
        return f'$\mathsf{{{-t:4.0f} ms before merger}}$'.replace(' ', '\ ')
    return f'$\mathsf{{{t:4.0f} ms after merger}}$'.replace(' ', '\ ')

def roundup(x, to):
    return int(math.ceil(x / to)) * to

def flux_astro(espl):
    I_0 = 1e-12*u.W/u.m**2
    D = 300 * u.Mpc
    flux = np.power(10, espl/10) * I_0 * D**2
    return flux.to("Lsun cm2 cm-2").value

directory = os.path.splitext(args.state)[0]
output_dir = directory
if args.output_dir is not None:
    output_dir = args.output_dir

# Create output directory if it doesn't exist
try:
    os.mkdir(output_dir)
except:
    pass


if args.m1 is None:
    raise "Primary mass is not specified"
if args.m2 is None:
    raise "Secondary mass is not specified"
if args.merger_time is None:
    raise "Merger time is not specified"

total_mass = args.m1 + args.m2

gu = GU(total_mass)

espl_limits = (20, 110)

first_frame = plt.imread(f'{directory}/frame.{0:06d}.png')
print(first_frame.shape)

px = 1/plt.rcParams['figure.dpi'] 
figsize = (first_frame.shape[1]*px, first_frame.shape[0]*px)

# disable automatic camera reset on 'Show'
pv._DisableFirstRenderCameraReset()

print(f'[{args.task_id:04d}] Loading state...')
# load state
pv.LoadState(args.state)
print(f'[{args.task_id:04d}] Loaded state.')

# get animation scene
animation = pv.GetAnimationScene()

# Calculate timestamps
global_frame_times = np.arange(
    start = animation.StartTime,
    stop = animation.EndTime,
    step = args.frame_spacing,
    )
animation.NumberOfFrames = len(global_frame_times)
split_global_frame_times = np.array_split(global_frame_times, args.total_task_number)

frame_number_offset = 0
for a in split_global_frame_times[:args.task_id]:
    frame_number_offset += len(a)

frame_times = split_global_frame_times[args.task_id]

# Set time
for i, frame_time in enumerate(frame_times):
    global_frame_id = frame_number_offset + i
    print(f'[{args.task_id:04d}] Postprocessing frame #{global_frame_id:06d} ({1+i:04d} out of local batch of size {len(frame_times):04d})')
    
    fig = plt.figure(figsize=figsize)
    ax = fig.gca()
    fig.tight_layout()

    frame = plt.imread(f'{directory}/frame.{global_frame_id:06d}.png')

    ax.imshow(frame)
    ax.margins(0, 0)
    ax.axis('off')

    ## Left upper corner
    # Masses
    ax.annotate(f"$\mathsf{{M = {args.m1:.0f}\ M_{{\odot}} + {args.m2:.0f}\ M_{{\odot}}}}$",
        xy=(0.02, 0.98), xycoords='axes fraction', fontsize=14,
        xytext=(5, -10), textcoords='offset points',
        color = 'white'
    )

    # Left bottom corner
    # Time
    ax.annotate(time_text(frame_time, mass=total_mass, merger=args.merger_time),
        xy=(0.02, 0.02), xycoords='axes fraction', fontsize=14,
        xytext=(5, 10), textcoords='offset points',
        color = 'white'
    )

    if args.astro_units_flux is False:
        # Plot eSPL
        norm = mpl.colors.Normalize(vmin=espl_limits[0], vmax=espl_limits[1])
        mappable = mpl.cm.ScalarMappable(norm=norm, cmap='viridis')
        cbaxes = ax.inset_axes(bounds=[0.94, 0.035, 0.01, 0.4])
        cb = fig.colorbar(mappable, cax=cbaxes, orientation='vertical')
        cbaxes.yaxis.label.set_color('black')
        cbaxes.tick_params(axis='y', colors='white', labelsize=12)
        cbaxes.yaxis.set_major_formatter(FuncFormatter(
            lambda x, pos: f'$\mathsf{{{x:.0f}}}$'
            ))
        for spine in cbaxes.spines:
            cbaxes.spines[spine].set_color('#3e3e3e')
        ax.annotate("$\mathsf{eSPL\ [dB]}$",
            xy=(0.915, 0.47), xycoords='axes fraction', fontsize=14,
            xytext=(5, 10), textcoords='offset points',
            color = 'white'
        )
        ax.annotate("$\mathsf{D_L=300\ Mpc}$",
            xy=(0.915, 0.44), xycoords='axes fraction', fontsize=9,
            xytext=(5, 10), textcoords='offset points',
            color = 'white'
        )
    else:
        # Plot L_sun/sr
        flux_astro_log_limits = (
            np.log10(flux_astro(espl_limits[0])),
            np.log10(flux_astro(espl_limits[1])),
        )
        norm = mpl.colors.Normalize(vmin=flux_astro_log_limits[0], vmax=flux_astro_log_limits[1])
        mappable = mpl.cm.ScalarMappable(norm=norm, cmap='viridis')
        cbaxes = ax.inset_axes(bounds=[0.94, 0.035, 0.01, 0.4])
        cb = fig.colorbar(mappable, cax=cbaxes, orientation='vertical',
                          ticks=np.arange(np.ceil(flux_astro_log_limits[0]), np.ceil(flux_astro_log_limits[1]), 1),
                          )
        cbaxes.yaxis.label.set_color('black')
        cbaxes.tick_params(axis='y', colors='white', labelsize=12)
        cbaxes.yaxis.set_major_formatter(FuncFormatter(
            lambda x, pos: f'$\mathsf{{ 10^{{ {x:.0f} }} }}$'
            ))
        for spine in cbaxes.spines:
            cbaxes.spines[spine].set_color('#3e3e3e')
        ax.annotate("$\mathsf{F\ [L_{\odot}\ sr^{-1}]}$",
            xy=(0.905, 0.44), xycoords='axes fraction', fontsize=14,
            xytext=(5, 10), textcoords='offset points',
            color = 'white'
        )

    output_filename = f'{directory}/frame-overplotted.{global_frame_id:06d}.png'
    if args.astro_units_flux is True:
        output_filename = f'{directory}/frame-overplotted-astro_units_flux.{global_frame_id:06d}.png'

    fig.savefig(output_filename, pad_inches=0, bbox_inches='tight')
    plt.close(fig)
