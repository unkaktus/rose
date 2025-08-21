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
    # energy_flux_source = pv.FindSource('EnergyFlux')
    # if energy_flux_source is None:
    #     raise "No source named EnergyFlux found"
    # domain_size = energy_flux_source.Size
    return 400
    # return domain_size

def roundz(x, n):
    return round(x, n)+0.0


## XXX Adapted
def time_text(time, mass, merger) -> str:
    t = mass * GU(1).T.to('s').value * (time - merger)
    if t < 0:
        return f'$\mathsf{{{t:4.3f}~s}}$'.replace(' ', '~')
    return ""
    # return f'$\mathsf{{+{t:4.3f}~s}}$'.replace(' ', '~')

def camera_angles():
    renderView = pv.GetActiveViewOrCreate('RenderView')
    camera_position = renderView.CameraPosition
    _, camera_theta, camera_phi = astropy.coordinates.cartesian_to_spherical(camera_position[0], camera_position[1], camera_position[2])
    camera_theta = np.pi/2 * u.rad  - camera_theta
    # camera_theta = 1.4 * u.rad
    return camera_theta.value, camera_phi.value

def find_value(timeseries, time):
    index = np.argmin(np.abs(time - timeseries[:,0]))
    return timeseries[index,1:]