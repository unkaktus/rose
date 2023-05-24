import h5py
import os
import numpy as np
from utils import get_data_root

def read_spins():
    chi_timeseries = []
    data_root = get_data_root()
    horizons_file = h5py.File(os.path.join(data_root, "Horizons.h5"))
    for key in horizons_file.keys():
        chi_timeseries.append(horizons_file[key]['chiInertial.dat'][:])
    horizons_file.close()
    return chi_timeseries

def spin_projections(chi):
    L = np.array([0, 0, 1])
    chi_parallel = np.dot(chi, L)
    chi_perp = np.linalg.norm(chi - chi_parallel*L)
    return chi_parallel, chi_perp