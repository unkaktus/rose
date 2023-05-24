import paraview.simple as pv
import os
import astropy
import numpy as np
from geounits import *

def get_data_root():
    energy_flux_source = pv.FindSource('EnergyFlux')
    if energy_flux_source is None:
        raise "No source named EnergyFlux found"
    waveform_file = f'{energy_flux_source.FileName}'
    return os.path.dirname(waveform_file)

def get_domain_size():
    energy_flux_source = pv.FindSource('EnergyFlux')
    if energy_flux_source is None:
        raise "No source named EnergyFlux found"
    domain_size = energy_flux_source.Size
    return domain_size

def roundz(x, n):
    return round(x, n)+0.0

def time_text(time, mass, merger) -> str:
    t = GU(mass).T.to('ms').value * (time - merger)
    if t < 0:
        return f'$\mathsf{{{t:4.0f} ms}}$'.replace(' ', '~')
    return f'$\mathsf{{+{t:4.0f} ms}}$'.replace(' ', '~')

def camera_angles():
    renderView = pv.GetActiveViewOrCreate('RenderView')
    camera_position = renderView.CameraPosition
    _, camera_theta, camera_phi = astropy.coordinates.cartesian_to_spherical(camera_position[0], camera_position[1], camera_position[2])
    return camera_theta, camera_phi

def find_value(timeseries, time):
    index = np.argmin(np.abs(time - timeseries[:,0]))
    return timeseries[index,1:]