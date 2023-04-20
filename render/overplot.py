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
from astropy import units as u, constants as const
import astropy.coordinates
import toml
from attrdict import AttrDict

import scri
import spherical
import quaternionic


berlin = {
    "U1": "#62AD2D",
    "U2": "#E94D10",
    "U3": "#00A192",
    "U4": "#FFD401",
    "U5": "#815237",
    "U6": "#846DAA",
    "U7": "#009AD9",
    "U8": "#005A99",
    "U9": "#F18800",
    "S9": "#9B2B48",
    "S7": "#846DAA",
    "S8": "#62AD2D",
    "S2": "#007B3D",
}

# better pictures and legends
plt.rc('figure', dpi=200)
plt.rc('text', usetex=True)
plt.style.use('dark_background')
plt.rc("savefig", facecolor="#121212")
plt.rc("savefig", edgecolor="#121212")


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

def evaluate_swsh(spin_weight, ell_max, coefficients, theta, phi):
    wigner = spherical.Wigner(ell_max=ell_max)
    direction = quaternionic.array.from_spherical_coordinates(theta, phi)
    total = np.zeros(coefficients.shape[0], dtype=complex)

    for l in range(abs(spin_weight), ell_max + 1):
        for abs_m in range(0, l + 1):
            for sign_m in (-1, 1):
                m = abs_m * sign_m
                mode_profile = wigner.sYlm(s=spin_weight, R=direction)[spherical.LM_index(l, m, 0)]
                total += coefficients[:,spherical.LM_index(l, m, 2)] * mode_profile
    return total

directory = os.path.splitext(args.state)[0]
output_dir = directory
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


# TODO: take these limits from the state file
espl_limits = (40, 110)

first_frame = plt.imread(f'{directory}/frame.{0:06d}.png')
print(first_frame.shape)

px = 1/plt.rcParams['figure.dpi'] 
# figsize = (first_frame.shape[1]*px, first_frame.shape[0]*px)
figsize = (first_frame.shape[1]*px, 1.2*first_frame.shape[0]*px)

# disable automatic camera reset on 'Show'
pv._DisableFirstRenderCameraReset()

print(f'[{args.task_id:04d}] Loading state...')
# load state
pv.LoadState(args.state)
print(f'[{args.task_id:04d}] Loaded state.')

# Get the waveform file
energy_flux_source = pv.FindSource('EnergyFlux')
if energy_flux_source is None:
    raise "No source named EnergyFlux found"
waveform_file = f'{energy_flux_source.FileName}/{energy_flux_source.Subfile}'
strain = scri.SpEC.read_from_h5(waveform_file)

merger_time = strain.t[np.argmax(np.abs(strain.data[:,spherical.LM_index(2, 2, 2)]))]
print(f'Merger time: {merger_time}')

# Get camera position
renderView = pv.GetActiveViewOrCreate('RenderView')
camera_position = renderView.CameraPosition
_, camera_theta, camera_phi = astropy.coordinates.cartesian_to_spherical(camera_position[0], camera_position[1], camera_position[2])

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
    gs = GridSpec(2, 1, figure=fig,
                  height_ratios=[5, 1],
                  width_ratios=[1],
                  hspace=0.0
                )
    fig.tight_layout()

    ax = fig.add_subplot(gs[0, :])

    frame = plt.imread(f'{directory}/frame.{global_frame_id:06d}.png')

    ax.imshow(frame)
    ax.margins(0, 0)
    ax.axis('off')

    ## Left upper corner
    # Masses
    ax.annotate(f"$\mathsf{{M = {datacorner.M:.0f}\ M_{{\odot}} }}$",
        xy=(0.02, 0.98), xycoords='axes fraction', fontsize=14,
        xytext=(5, -10), textcoords='offset points',
        color = 'white'
    )
    ax.annotate(f"$\mathsf{{q = 1:{datacorner.q:.0f} }}$",
        xy=(0.02, 0.95), xycoords='axes fraction', fontsize=14,
        xytext=(5, -10), textcoords='offset points',
        color = 'white'
    )
    ax.annotate(f"$\mathsf{{ \\vec{{\chi_1}} = ({datacorner.chi1[0]:.1f}, {datacorner.chi1[1]:.1f}, {datacorner.chi1[2]:.1f}) }}$",
        xy=(0.02, 0.91), xycoords='axes fraction', fontsize=14,
        xytext=(5, -10), textcoords='offset points',
        color = 'white'
    )
    ax.annotate(f"$\mathsf{{ e = {datacorner.e:.1f} }}$",
        xy=(0.02, 0.88), xycoords='axes fraction', fontsize=14,
        xytext=(5, -10), textcoords='offset points',
        color = 'white'
    )


    # Left bottom corner
    # Time
    ax.annotate(time_text(frame_time, mass=datacorner.M, merger=merger_time),
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

    # Now add the waveform
    h = evaluate_swsh(spin_weight=-2, ell_max=8, coefficients=strain.data, theta=camera_theta, phi=camera_phi)
    waveform_ax = fig.add_subplot(gs[1, :])
    h_plus, h_cross = np.real(h), np.imag(h)
    waveform_ax.plot(strain.t, h_cross, color="gray")
    waveform_ax.plot(strain.t, h_plus, color="darkgray")

    mask = strain.t < frame_time
    waveform_ax.plot(strain.t[mask], h_cross[mask], color=berlin['S2'], label=r'$\mathsf{h_{\times}}$')
    waveform_ax.plot(strain.t[mask], h_plus[mask], color=berlin['S8'], label=r'$\mathsf{h_{+}}$')
    waveform_ax.axvline(x=frame_time, color="white")
    waveform_ax.axis('off')
    waveform_ax.set_xlim(0,strain.t[-1])
    h_max = np.max(np.abs(h))
    waveform_ax.set_ylim(-1.2*h_max, 1.2*h_max)
    waveform_ax.legend(fontsize=14, loc='upper left', ncol=3, frameon=False, labelspacing=0, handletextpad=0.4, borderpad=0.2, borderaxespad=1.5)

    output_filename = f'{output_dir}/frame-overplotted.{global_frame_id:06d}.png'
    if args.astro_units_flux is True:
        output_filename = f'{output_dir}/frame-overplotted-astro_units_flux.{global_frame_id:06d}.png'

    fig.savefig(output_filename, pad_inches=0, bbox_inches='tight')
    plt.close(fig)

print(f'[{args.task_id:04d}] Done.')