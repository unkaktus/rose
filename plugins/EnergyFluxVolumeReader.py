
import os
import hashlib
import logging
import tempfile
import time
import concurrent.futures

import numpy as np
import spherical
import quaternionic
from spherical_functions import LM_index
import scri
from astropy import units as u, constants as const
import scipy
import psutil


from vtkmodules.vtkCommonDataModel import vtkUniformGrid
from vtkmodules.vtkCommonCore import vtkDataArraySelection
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.numpy_interface import dataset_adapter as dsa
from paraview.vtk.util import numpy_support as vtknp
from paraview.util.vtkAlgorithm import smdomain, smhint, smproperty, smproxy
from paraview import util


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

rose_cache_dir = os.environ.get("ROSE_CACHE_DIR")
if rose_cache_dir is None:
    rose_cache_dir = tempfile.mkdtemp(prefix='rose-')

workers = int(os.environ.get("ROSE_WORKERS"))
if workers is None:
    workers = psutil.cpu_count(logical=False)

EnergyFluxLog10ArrayName = "Energy flux (log10)"
ESPLArrayName = "eSPL in dB"

# Needed by ParaView
def create_modified_callback(anobject):
    import weakref

    weakref_obj = weakref.ref(anobject)
    anobject = None

    def _markmodified(*args, **kwars):
        o = weakref_obj()
        if o is not None:
            o.Modified()

    return _markmodified

def set_output_array(output, name, array):
    quantity_vtk = vtknp.numpy_to_vtk(array, deep=True)
    quantity_vtk.SetName(name)
    output.GetPointData().AddArray(quantity_vtk)


def energy_flux_hash(waveform_filename):
    key = f'energy_flux|{waveform_filename}'.encode('utf-8')
    return hashlib.sha256(key).hexdigest()[:16]


class GU():
    def __init__(self, v):
        self.v = v*u.Msun
        self.M = self._to(u.g)
        self.L = self._to(u.cm)
        self.T = self._to(u.s)
        self.Energy = self.MLT(1, 2, -2)
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

def eSPL(energy_flux, distance):
    gu = GU(1)
    I_0 = 1e-12*u.W/u.m**2
    intensity = energy_flux * gu.Luminosity / (distance**2)
    SPL = 10 * np.log10( (intensity/I_0).to('').value)
    return SPL


def get_timestep(algorithm, logger=None):
    if logger is None:
        logger = logging.getLogger(__name__)
    logger.debug("Getting current timestep...")
    executive = algorithm.GetExecutive()
    outInfo = executive.GetOutputInformation(0)
    if not outInfo.Has(executive.UPDATE_TIME_STEP()):
        logger.debug("No `UPDATE_TIME_STEP` found. Return 0.")
        return 0.0
    timestep = outInfo.Get(executive.UPDATE_TIME_STEP())
    logger.debug(f"Found `UPDATE_TIME_STEP`: {timestep}")
    return timestep


def set_timesteps(algorithm, timesteps, logger=None):
    if logger is None:
        logger = logging.getLogger(__name__)
    executive = algorithm.GetExecutive()
    outInfo = executive.GetOutputInformation(0)
    outInfo.Remove(executive.TIME_STEPS())
    outInfo.Remove(executive.TIME_RANGE())

    if timesteps is not None:
        for timestep in timesteps:
            outInfo.Append(executive.TIME_STEPS(), timestep)
        outInfo.Append(executive.TIME_RANGE(), timesteps[0])
        outInfo.Append(executive.TIME_RANGE(), timesteps[-1])
    logger.debug(
        f"Set data timesteps to {outInfo.Get(executive.TIME_RANGE())}."
    )


def get_mode_name(l, abs_m):
    return "({}, {}) Mode".format(l, abs_m)


class SWSHParameters:
    size: float
    num_points: float
    spin_weight: int
    ell_max: int

class SphericalGrid:
    r: np.ndarray
    theta: np.ndarray
    phi: np.ndarray


class CartesianGrid:
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray

    def __init__(self, p: SWSHParameters) -> None:
        X = np.linspace(-p.size, p.size, p.num_points)
        Y = X
        Z = X
        self.x, self.y, self.z = map(
            lambda arr: arr.flatten(order="F"), np.meshgrid(
                X, Y, Z, indexing="ij")
        )

    def spherical(self) -> SphericalGrid:
        spherical_grid = SphericalGrid()
        spherical_grid.r = np.sqrt(self.x**2 + self.y**2 + self.z**2)
        spherical_grid.theta = np.arccos(self.z / spherical_grid.r)
        spherical_grid.phi = np.arctan2(self.y, self.x)
        return spherical_grid

def smoothcos(x, v, scale, width):
    mask = x > scale + width/2
    v[mask] = 0
    mask = np.logical_and(scale - width/2 < x, x < scale + width/2)
    v[mask] *= np.cos((np.pi/2)* (x[mask]-scale+width/2) / width )**2
    return v


def Ylm(l, m, theta, phi, out=None):
    return scipy.special.sph_harm(m, l, phi, theta, out)


def all_modes(r, theta, phi, energy_flux, phase, timesteps, ell_max):
    # Compute quantity in the volume from the input waveform data
    quantity = np.zeros_like(r, dtype=complex)
    mode = np.zeros_like(quantity, dtype=complex)
    c_lm = np.zeros_like(quantity, dtype=complex)
    for l in range(0, ell_max + 1):
        for m in range(-l, l+1):
            Ylm(l=l, m=m, phi=phi, theta=theta, out=mode)
            c_lm_timeseries = np.array(energy_flux[:, LM_index(l,m,0)])

            # Interpolate c_lm using linear interpolation
            c_lm = np.interp(
                phase,
                timesteps,
                c_lm_timeseries,
                left=0.0,
                right=0.0,
            )
            mode *= c_lm
            quantity += mode
    del mode
    del c_lm
    return quantity


@smproxy.reader(
    name="EnergyFluxVolumeReader",
    label="Energy Flux Volume Reader",
    extensions="h5",
    file_description="HDF5 files",
)
class EnergyFluxVolumeReader(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self,
            nInputPorts=0,
            nOutputPorts=1,
            outputType="vtkUniformGrid",
        )

        self.energy_flux_cache_dir = os.path.join(rose_cache_dir, "energy_flux_cache")

        self._filename = None
        self._subfile = None

        self.spin_weight = 0
        self.ell_max = 8

        self.timesteps = None

        self.modes_selection = vtkDataArraySelection()
        for l in range(np.abs(self.spin_weight), self.ell_max + 1):
            for m in range(0, l + 1):
                self.modes_selection.AddArray(get_mode_name(l, m))
        self.modes_selection.AddObserver(
            "ModifiedEvent", create_modified_callback(self)
        )

        self.component_selection = vtkDataArraySelection()
        self.component_selection.AddArray(EnergyFluxLog10ArrayName)
        self.component_selection.AddArray(ESPLArrayName)
        self.component_selection.DisableAllArrays()
        self.component_selection.AddObserver(
            "ModifiedEvent", create_modified_callback(self)
        )

        self.component_selection.EnableArray(ESPLArrayName)


    @smproperty.stringvector(name="FileName")
    @smdomain.filelist()
    @smhint.filechooser(extensions="h5", file_description="HDF5 files")
    def SetFileName(self, value):
        self._filename = value
        self.Modified()

    @smproperty.stringvector(
        name="Subfile", default_values=["Extrapolated_N2.dir"]
    )

    def SetSubfile(self, value):
        self._subfile = value
        self.Modified()


    def _get_waveform_data(self):
        return self.waveform_data

    @smproperty.dataarrayselection(name="Modes")
    def GetModes(self):
        return self.modes_selection

    @smproperty.dataarrayselection(name="Components")
    def GetPolarizations(self):
        return self.component_selection

    @smproperty.doublevector(name="Size", default_values=100)
    def SetSize(self, value):
        self.size = value
        self.Modified()

    @smproperty.intvector(name="SpatialResolution", default_values=100)
    def SetSpatialResolution(self, value):
        self.num_points_per_dim = value
        self.Modified()


    @smproperty.intvector(name="EllMax", default_values=8)
    def SetEllMax(self, value):
        self.ell_max = value
        self.Modified()

    @smproperty.doublevector(name="TimeShift", default_values=0.0)
    def SetTimeShift(self, value):
        self.time_shift = value
        self.Modified()

    @smproperty.doublevector(name="ValueThreshold", default_values=1e-16)
    def SetValueThreshold(self, value):
        self.value_threshold = value
        self.Modified()

    @smproperty.doublevector(name="DistanceInMpc", default_values=300)
    def SetDistance(self, value):
        self.distance = value
        self.Modified()

    def _get_timesteps(self):
        ts = self.timesteps
        return np.linspace(ts[0], ts[-1], len(ts))

    @smproperty.doublevector(
        name="TimestepValues",
        information_only="1",
        si_class="vtkSITimeStepsProperty",
    )

    def _get_array_selection(self):
        return self._arrayselection_cd, self._arrayselection_pd

    def GetTimestepValues(self):
        timesteps = self._get_timesteps()
        set_timesteps(self, self._get_timesteps(), logger=logger)
        return timesteps.tolist() if timesteps is not None else None

    def RequestInformation(self, request, inInfo, outInfo):
        logger.debug("Requesting information...")

        self.load_data()

        info = outInfo.GetInformationObject(0)

        N = self.num_points_per_dim
        N_y = N
        N_z = N
        grid_extents = [0, N - 1, 0, N_y - 1, 0, N_z - 1]
        util.SetOutputWholeExtent(self, grid_extents)

        # This needs the time data from the waveform file, so we may have to
        # set the `TIME_RANGE` and `TIME_STEPS` already in the
        # WaveformDataReader.
        set_timesteps(self, self._get_timesteps(), logger=logger)

        D = self.size

        grid_params = SWSHParameters()
        grid_params.size = D
        grid_params.num_points = N
        grid_params.spin_weight = self.spin_weight
        grid_params.ell_max = self.ell_max

        start = time.time()
        grid = CartesianGrid(grid_params)
        self.spherical_grid = grid.spherical()
        end = time.time()
        logger.debug(f"Done grid creation in {end - start}")

        return 1

    def load_data(self):
        if self.timesteps is not None:
            return
        if self._filename is None:
            return
        if self._subfile is None:
            return
        
        abd = scri.SpEC.create_abd_from_h5("SXS", h = f"{self._filename}/{self._subfile}")
            
        # Get the timesteps
        self.timesteps = abd.t

        # Cache calculated energy flux
        energy_flux = np.array([])
        energy_flux_filename = os.path.join(
            self.energy_flux_cache_dir,
            f'{energy_flux_hash(self._filename)}.npy'
            )

        try:
            energy_flux = np.load(energy_flux_filename)
        except:
            # Calculate the energy flux modes coefficients
            energy_flux = (1/(4*np.pi)) * np.array(abd.sigma.dot * abd.sigma.dot.bar)
            if not os.path.exists(self.energy_flux_cache_dir):
                os.makedirs(self.energy_flux_cache_dir)
            np.save(energy_flux_filename, energy_flux)

        self.energy_flux = energy_flux

    def RequestData(self, request, inInfo, outInfo):
        logger.info("Requesting data...")

        set_timesteps(self, self._get_timesteps(), logger=logger)

        output = dsa.WrapDataObject(vtkUniformGrid.GetData(outInfo))

        t = get_timestep(self, logger=logger)
        N = self.num_points_per_dim
        D = self.size

        dx = 2.0 * D / (N-1)
        N_y = N
        N_z = N
        output.SetDimensions(N, N_y, N_z)
        output.SetOrigin(-D, -D, -D)
        output.SetSpacing(dx, dx, dx)


        start = time.time()
        # Compute scaled waveform phase on the grid
        phase = (t + self.time_shift) - self.spherical_grid.r
        end = time.time()
        logger.debug(f"Done phase calculation in {end - start}")

        start = time.time()
        # Split to workers
        logger.debug(f'workers: {workers}')
        r_w = np.array_split(self.spherical_grid.r, workers)
        theta_w = np.array_split(self.spherical_grid.theta, workers)
        phi_w = np.array_split(self.spherical_grid.phi, workers)
        phase_w = np.array_split(phase, workers)
        quantity_w = np.zeros(shape=(workers,), dtype=np.ndarray)

        def all_modes_w(worker):
            quantity_w_local = all_modes(r=r_w[worker],
                    theta=theta_w[worker],
                    phi=phi_w[worker],
                    energy_flux=self.energy_flux,
                    phase=phase_w[worker],
                    timesteps=self.timesteps,
                    ell_max=self.ell_max)
            logger.debug(f"{worker} Calculated all modes")
            return quantity_w_local

        logger.debug(f"creating {workers} workers")
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(all_modes_w, worker): worker for worker in range(workers)}
            for future in concurrent.futures.as_completed(futures):
                worker = futures[future]
                logger.debug(f'worker {worker} is done')
                try:
                    quantity_w[worker] = future.result()
                except Exception as exc:
                    print('%r generated an exception: %s' % (futures[future], exc))

        logger.debug(f"Executor is done")
        quantity = np.concatenate(quantity_w)

        del quantity_w, r_w, theta_w, phi_w, phase_w

        end = time.time()
        logger.debug(f"Done calculating modes in {end - start}")

        energy_flux = np.real(quantity)

        # Clip from below
        np.maximum(energy_flux, self.value_threshold, out=energy_flux)

        # Add entire sum to the output
        if self.component_selection.ArrayIsEnabled(EnergyFluxLog10ArrayName):
            energy_flux_log = np.log10(energy_flux)
            set_output_array(output, EnergyFluxLog10ArrayName, energy_flux_log)

        # Effective Sound Pressure Level in dB
        if self.component_selection.ArrayIsEnabled(ESPLArrayName):
            espl = eSPL(energy_flux, distance=self.distance*u.Mpc)
            espl = smoothcos(self.spherical_grid.r, espl, scale=self.size*0.95, width=5)
            set_output_array(output, ESPLArrayName, espl)


        return 1
