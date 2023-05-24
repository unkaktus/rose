import spherical
import quaternionic
import h5py
import scri
import numpy as np
import paraview.simple as pv


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

def read_strain():
    energy_flux_source = pv.FindSource('EnergyFlux')
    if energy_flux_source is None:
        raise "No source named EnergyFlux found"
    waveform_file = f'{energy_flux_source.FileName}/{energy_flux_source.Subfile}'
    strain = scri.SpEC.read_from_h5(waveform_file)
    return strain

def strain_merger_time(h):
    merger_time = h.t[np.argmax(np.abs(h.data[:,spherical.LM_index(2, 2, 2)]))]
    return merger_time
