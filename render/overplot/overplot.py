#!/usr/bin/env -S python -u

import argparse

parser = argparse.ArgumentParser('overplot.py')
parser.add_argument("--total-task-number", type=int, help="Total number of tasks", default=1)
parser.add_argument("--task-id", type=int, help="Current task ID", default=0)
parser.add_argument("--state", type=str, help="State filename")
parser.add_argument("--output-dir", type=str, help="Path to output directory", default=None)

parser.add_argument("--frame-spacing", type=float, help="Frame spacing in total masses", default=10)
parser.add_argument("--astro-units-flux", type=bool, help="Use astrophysical units for flux (L_sun/sr)", default=False)

args = parser.parse_args()

import os
import math
import numpy as np
import paraview.simple as pv

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FuncFormatter
import astropy.coordinates
import toml
from attrdict import AttrDict

from geounits import *
from colorbar import *
from colors import *
from spins import *
from waveform import *
from utils import *

# better pictures and legends
plt.rc('figure', dpi=200)
plt.rc('text', usetex=True)
plt.style.use('dark_background')
plt.rc("savefig", facecolor=background_color)
plt.rc("savefig", edgecolor=background_color)

aspect_ratio = 16/9
fontsize = 20

aei_logo_filename = "/home/imarkin/rose/logos/aei.png"
sxs_logo_filename = "/home/imarkin/rose/logos/sxs.png"

directory = os.path.splitext(args.state)[0]
source_dir = os.path.join(directory, 'rendered')
zoom_dir = os.path.join(f'{directory}_zoom', 'rendered')

output_dir = os.path.join(directory, 'overplotted')

if args.output_dir is not None:
    output_dir = args.output_dir

# Create output directory if it doesn't exist
try:
    os.mkdir(output_dir)
except:
    pass

# Load datacorner.toml file
datacorner = AttrDict(toml.load(
    os.path.join(os.path.dirname(directory), "datacorner.toml")
))

if datacorner.M is None:
    raise "Total mass M is not specified"
if datacorner.q is None:
    raise "Mass ratio q is not specified"

first_frame = plt.imread(f'{source_dir}/frame.{0:06d}.png')

px = 1/plt.rcParams['figure.dpi'] 
figheight = 1.2 * first_frame.shape[0]*px
figsize = (aspect_ratio*figheight, figheight)

# disable automatic camera reset on 'Show'
pv._DisableFirstRenderCameraReset()

print(f'[{args.task_id:04d}] Loading state...')
# load state
pv.LoadState(args.state)
print(f'[{args.task_id:04d}] Loaded state.')

strain = read_strain()

merger_time = strain_merger_time(strain)
print(f'Merger time: {merger_time}')

# Get camera position
camera_theta, camera_phi = camera_angles()
# get animation scene
animation = pv.GetAnimationScene()

# Calculate timestamps
input_frame_times = np.arange(
    start = animation.StartTime,
    stop = animation.EndTime,
    step = args.frame_spacing,
    )
# Save the IDs of the frames
input_frame_ids = np.arange(len(input_frame_times))

animation.NumberOfFrames = len(input_frame_times)


global_frame_times = input_frame_times
global_frame_ids = input_frame_ids

# Always let the waves propagate to the domain edge
start_time = get_domain_size()

# Skip number of frames
if datacorner.start_time is not None:
    start_time += datacorner.start_time
    # print(f'Skipping to t={datacorner.start_time}')
    mask = input_frame_times > start_time
    global_frame_times = input_frame_times[mask]
    global_frame_ids = input_frame_ids[mask]

    mask = strain.t > start_time
    strain.t = strain.t[mask]
    strain.data = strain.data[mask]



# Read eccentricity
e_timeseries = np.array([[0, datacorner.e]]) # Fill with constant eccentricity by default
if hasattr(datacorner, 'e_file'):
    e_timeseries = np.loadtxt(datacorner.e_file)

split_global_frame_times = np.array_split(global_frame_times, args.total_task_number)
split_global_frame_ids = np.array_split(global_frame_ids, args.total_task_number)

frame_number_offset = 0
for a in split_global_frame_times[:args.task_id]:
    frame_number_offset += len(a)

# Get task-local frames
frame_times = split_global_frame_times[args.task_id]
frame_ids = split_global_frame_ids[args.task_id]

def add_border_to_axes(ax, color='darkgrey'):
    ax.margins(0, 0)
    ax.tick_params(
        axis='both',
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,
        left=False,         # ticks along the top edge are off
        labelbottom=False,
        labelleft=False,
        ) # labels along the bottom edge are off
    for spine in ax.spines:
        ax.spines[spine].set_color(color)

def axis_ratio(ax):
    return ax.bbox.width / ax.bbox.height

def imshow(ax, filename, zoom=1):
    image = plt.imread(filename)
    ax.imshow(image, origin='upper')
    if zoom != 1:
        (center_x, center_y) = (image.shape[1]//2, image.shape[0]//2)
        extent = (image.shape[1]//zoom, image.shape[0]//zoom)
        ax.set_xlim((center_x - axis_ratio(ax)*extent[0]//2, center_x + axis_ratio(ax)*extent[0]//2))
        ax.set_ylim((center_y - extent[1]//2, center_y + extent[1]//2))

def add_render(fig, axes, source_dir, frame_id):
    ax = fig.add_subplot(axes)
    filename = f'{source_dir}/frame.{frame_id:06d}.png'
    imshow(ax, filename, zoom=1.5)
    ax.margins(0, 0)
    ax.axis('off')
    # add_border_to_axes(ax)

def add_zoom_in(fig, axes, frame_id):
    ax = fig.add_subplot(axes)
    ax.margins(0, 0)
    add_border_to_axes(ax)
    print(axis_ratio(ax))
    filename = f'{zoom_dir}/frame.{frame_id:06d}.png'
    imshow(ax, filename)

def add_logos(fig, logos):
    for logo in logos:
        (axes, filename) = logo
        ax = fig.add_subplot(axes)
        ax.margins(0, 0)
        image = plt.imread(filename)
        ax.imshow(image)
        ax.axis('off')
        # add_border_to_axes(ax)

def add_colorbar(fig, axes):
    ax = fig.add_subplot(axes)
    ax.margins(0, 0)
    ax.axis('off')
    # add_border_to_axes(ax)
    cbar_bounds = [0.22, 0.05, 0.03, 0.9]

    if args.astro_units_flux is False:
        # Plot eSPL
        norm = mpl.colors.Normalize(vmin=datacorner.espl_limits[0], vmax=datacorner.espl_limits[1])
        mappable = mpl.cm.ScalarMappable(norm=norm, cmap='viridis')
        cbaxes = ax.inset_axes(bounds=cbar_bounds)
        fig.colorbar(mappable, cax=cbaxes, orientation='vertical')
        cbaxes.yaxis.label.set_color('black')
        cbaxes.tick_params(axis='y', colors='white', labelsize=12)
        cbaxes.yaxis.set_major_formatter(FuncFormatter(
            lambda x: f'$\mathsf{{{x:.0f}}}$'
            ))
        for spine in cbaxes.spines:
            cbaxes.spines[spine].set_color(background_color)
        ax.annotate("$\mathsf{eSPL\ [dB]}$",
            xy=(0.915, 0.47), xycoords='axes fraction', fontsize=fontsize,
            color = 'white'
        )
        ax.annotate("$\mathsf{D_L=300\ Mpc}$",
            xy=(0.915, 0.44), xycoords='axes fraction', fontsize=fontsize*2/3,
            xytext=(5, 10), textcoords='offset points',
            color = 'white'
        )
    if args.astro_units_flux is True:
        # Plot L_sun/sr
        flux_astro_log_limits = (
            np.log10(flux_astro(datacorner.espl_limits[0])),
            np.log10(flux_astro(datacorner.espl_limits[1])),
        )
        norm = mpl.colors.Normalize(vmin=flux_astro_log_limits[0], vmax=flux_astro_log_limits[1])
        mappable = mpl.cm.ScalarMappable(norm=norm, cmap='viridis')
        cbaxes = ax.inset_axes(bounds=cbar_bounds)
        cb = fig.colorbar(mappable, cax=cbaxes, orientation='vertical',
                          ticks=np.arange(np.ceil(flux_astro_log_limits[0]), np.ceil(flux_astro_log_limits[1]), 1),
                          )
        cbaxes.yaxis.label.set_color('black')
        cbaxes.tick_params(axis='y', colors='white', labelsize=fontsize*3/4)
        cbaxes.yaxis.set_major_formatter(FuncFormatter(
            lambda x, pos: f'$\mathsf{{ 10^{{ {x:.0f} }} }}$'
            ))
        for spine in cbaxes.spines:
            cbaxes.spines[spine].set_color(background_color)
        ax.annotate("$\mathsf{F\ [L_{\odot}\ sr^{-1}]}$",
            xy=(cbar_bounds[0]-0.075, cbar_bounds[1] + cbar_bounds[3]+0.03), xycoords='axes fraction', fontsize=fontsize,
            # xytext=(5, 10), textcoords='offset points',
            color = 'white'
        )

def add_datacorner(fig, axes, frame_time):
    luc_ax = fig.add_subplot(axes)
    luc_ax.margins(0, 0)
    luc_ax.axis('off')
    x_position = 0.1
    y_position = 0.9
    line_height = 0.1
    xy = [x_position, y_position]
    # Masses
    luc_ax.annotate(f"$\mathsf{{M = {datacorner.M:.0f}\ M_{{\odot}} }}$",
        xy=xy, xycoords='axes fraction', fontsize=fontsize,
        color = 'white'
    )
    xy[1] -= line_height
    luc_ax.annotate(f"$\mathsf{{q = 1:{datacorner.q:.0f} }}$",
        xy=xy, xycoords='axes fraction', fontsize=fontsize,
        color = 'white'
    )

    chi_timeseries = read_spins()

    xy[1] -= line_height*1.2
    chi1 = find_value(chi_timeseries[0], frame_time)
    luc_ax.annotate(f"$\mathsf{{ \\vec{{\chi_1}} = ( {roundz(chi1[0], 1):.1f},~{roundz(chi1[1], 1):.1f},~{roundz(chi1[2], 1):.1f} ) }}$",
        xy=xy, xycoords='axes fraction', fontsize=fontsize,
        color = 'white'
    )
    xy[1] -= line_height*1.2
    chi2 = find_value(chi_timeseries[1], frame_time)
    luc_ax.annotate(f"$\mathsf{{ \\vec{{\chi_2}} = ( {roundz(chi2[0], 1):.1f},~{roundz(chi2[1], 1):.1f},~{roundz(chi2[2], 1):.1f} ) }}$",
        xy=xy, xycoords='axes fraction', fontsize=fontsize,
        color = 'white'
    )
    xy[1] -= line_height
    e = find_value(e_timeseries, frame_time-merger_time)[0]
    luc_ax.annotate(f"$\mathsf{{ e = {e:.1f} }}$",
        xy=xy, xycoords='axes fraction', fontsize=fontsize,
        color = 'white'
    )

def add_waveform(fig, axes, frame_time):
    waveform_ax = fig.add_subplot(axes)
    h = evaluate_swsh(spin_weight=-2, ell_max=8, coefficients=strain.data, theta=camera_theta, phi=camera_phi)
    h_max = np.max(np.abs(h))
    h_plus, h_cross = np.real(h), np.imag(h)
    limit_factor = 1.2
    
    waveform_ax.set_xlim(start_time,strain.t[-1])
    waveform_ax.set_ylim(-limit_factor*h_max, limit_factor*h_max)

    # Plot the gray future waveform
    waveform_ax.plot(strain.t, h_cross, color="gray")
    waveform_ax.plot(strain.t, h_plus, color="darkgray")

    # Plot the past waveform with colors
    mask = strain.t < frame_time
    waveform_ax.plot(strain.t[mask], h_cross[mask], color=berlin['S2'], label=r'$\mathsf{h_{\times}}$')
    waveform_ax.plot(strain.t[mask], h_plus[mask], color=berlin['S8'], label=r'$\mathsf{h_{+}}$')
    waveform_ax.axvline(x=frame_time,
                        ymax = 0.5 + 0.5*(2/3)/limit_factor,
                        color="white")
    waveform_ax.margins(0, 0)
    waveform_ax.axis('off')
    waveform_ax.legend(fontsize=fontsize, loc='upper left', ncol=3, frameon=False, labelspacing=0, handletextpad=0.4, borderpad=0.1, borderaxespad=0.1)

def add_clock(fig, axes, frame_time):
    lbc_ax = fig.add_subplot(axes)
    lbc_ax.margins(0, 0)
    lbc_ax.axis('off')
    # Time
    lbc_ax.annotate(time_text(frame_time, mass=datacorner.M, merger=merger_time),
        xy=(0.1, 0.1), xycoords='axes fraction', fontsize=fontsize,
        xytext=(5, 10), textcoords='offset points',
        color = 'white'
    )

# Main frame loop
for i, frame_time in enumerate(frame_times):
    global_frame_id = frame_number_offset + i
    print(f'[{args.task_id:04d}] Postprocessing frame #{global_frame_id:06d} ({1+i:04d} out of local batch of size {len(frame_times):04d})')
    
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(9*2, 16*2, figure=fig,
                  height_ratios=np.ones(9*2),
                  width_ratios=np.ones(16*2),
                  hspace=0.0,
                  wspace=0.0
                )
    axes = {
        "render": gs[0:16,  4:28],
        "zoomin": gs[0:6, 24:],
        "datacorner": gs[0:8, 0:8],
        "colorbar": gs[6:14, 0:6],
        "clock": gs[14:16, 0:8],
        "waveform": gs[16:18, 0:30],
        "logo1": gs[16:18, 30:31],
        "logo2": gs[16:18, 31:],
    }
    # gs has [y, x] indexing
    fig.tight_layout()

    add_render(fig, axes["render"], source_dir, frame_id=frame_ids[i])
    add_datacorner(fig, axes["datacorner"], frame_time)
    add_colorbar(fig, axes["colorbar"])
    add_clock(fig, axes["clock"], frame_time)
    add_waveform(fig, axes["waveform"], frame_time=frame_time)
    add_zoom_in(fig, axes["zoomin"], frame_id=frame_ids[i])
    add_logos(fig, [
        (axes["logo1"], aei_logo_filename),
        (axes["logo2"], sxs_logo_filename),
        ])
    
    gs.update(hspace=0.0, wspace=0.0)

    output_filename = f'{output_dir}/frame-overplotted.{global_frame_id:06d}.png'

    fig.savefig(output_filename, pad_inches=0, bbox_inches='tight')
    plt.close(fig)

print(f'[{args.task_id:04d}] Done.')