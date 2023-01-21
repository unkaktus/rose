import spherical
import quaternionic
from spherical_functions import LM_index
import scri
import numpy as np
import os
import hashlib
import time
import logging
import sys
import tempfile

from vtkmodules.vtkCommonDataModel import vtkUniformGrid
from vtkmodules.vtkCommonCore import vtkDataArraySelection
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.numpy_interface import dataset_adapter as dsa
from paraview.vtk.util import numpy_support as vtknp
from paraview.util.vtkAlgorithm import smdomain, smproperty, smproxy
from paraview import util


logger = logging.getLogger(__name__)
rose_cache_dir = os.environ.get("ROSE_CACHE_DIR", "")
if rose_cache_dir == "":
    rose_cache_dir = tempfile.mkdtemp(prefix='rose-')


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


def swsh_grid_hash(p: SWSHParameters):
    key = f'{p.size}|{p.num_points}|{p.spin_weight}|{p.ell_max}'.encode(
        'utf-8')
    return hashlib.sha256(key).hexdigest()[:16]


def load_swsh_grid(params: SWSHParameters, cache_dir: str):
    if cache_dir == "":
        return None
    grid_hash = swsh_grid_hash(params)
    cache_filename = os.path.join(
        cache_dir,
        f'{grid_hash}.npy',
    )
    if not os.path.exists(cache_filename):
        logger.info("Cache miss: {swsh_grid_cache_filename}")
        return None

    logger.info(
        f"Loading SWSH grid from file '{cache_filename}'..."
    )
    swsh_grid = np.load(cache_filename)
    return swsh_grid


def create_swsh_grid(p: SWSHParameters):
    logger.info("Creating SWSH grid...")

    cartesian_grid = CartesianGrid(p)
    spherical_grid = cartesian_grid.spherical()

    angles = quaternionic.array.from_spherical_coordinates(
        spherical_grid.theta, spherical_grid.phi)
    swsh_grid = spherical.Wigner(p.ell_max).sYlm(s=p.spin_weight, R=angles)
    return swsh_grid


def save_swsh_grid(swsh_grid, p: SWSHParameters, cache_dir: str):
    if cache_dir == "":
        return
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)
    grid_hash = swsh_grid_hash(p)
    cache_filename = os.path.join(
        cache_dir,
        f'{grid_hash}.npy',
    )
    if not os.path.exists(cache_filename):
        np.save(cache_filename, swsh_grid)
        logger.debug(f'SWSH grid cache saved to {cache_filename}.')


def load_or_create_swsh_grid(p: SWSHParameters, cache_dir: str):
    swsh_grid = load_swsh_grid(p, cache_dir)
    if swsh_grid is None:
        swsh_grid = create_swsh_grid(p)
        save_swsh_grid(swsh_grid, p, cache_dir)
    return swsh_grid

class GridAdjustParameters():
    scale_with_distance: bool


def adjust_swsh_grid(swsh_grid, grid_params: SWSHParameters, params: GridAdjustParameters):
    spherical_grid = CartesianGrid(grid_params).spherical()
    if params.scale_with_distance:
        swsh_grid /= (spherical_grid.r**2 + 1.0e-30).reshape(spherical_grid.r.shape + (1,))
    return swsh_grid, spherical_grid



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


@smproxy.filter(label="Energy Flux To Volume")
@smproperty.input(name="WaveformData", port_index=0)
@smdomain.datatype(dataTypes=["vtkTable"])
class EnergyFluxToVolume(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self,
            nInputPorts=1,
            nOutputPorts=1,
            outputType="vtkUniformGrid",
        )

        self.spin_weight = 0
        self.ell_max = 8

        self.modes_selection = vtkDataArraySelection()
        for l in range(np.abs(self.spin_weight), self.ell_max + 1):
            for m in range(0, l + 1):
                self.modes_selection.AddArray(get_mode_name(l, m))
        self.modes_selection.AddObserver(
            "ModifiedEvent", create_modified_callback(self)
        )

        self.component_selection = vtkDataArraySelection()
        self.component_selection.AddArray("Energy flux")
        self.component_selection.AddObserver(
            "ModifiedEvent", create_modified_callback(self)
        )

    def FillInputPortInformation(self, port, info):
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkTable")

    def _get_waveform_data(self):
        return dsa.WrapDataObject(self.GetInputDataObject(0, 0))

    @smproperty.dataarrayselection(name="Modes")
    def GetModes(self):
        return self.modes_selection

    @smproperty.intvector(name="StoreIndividualModes", default_values=False)
    def SetStoreIndividualModes(self, value):
        self.store_individual_modes = value
        self.Modified()

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

    @smproperty.intvector(name="ScaleWithDistance", default_values=False)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetScaleWithDistance(self, value):
        self.scale_with_distance = value
        self.Modified()

    @smproperty.doublevector(name="ValueThreshold", default_values=1e-16)
    def SetValueThreshold(self, value):
        self.value_threshold = value
        self.Modified()

    @smproperty.stringvector(name="SwshCacheDirectory", default_values=os.path.join(rose_cache_dir, "swsh_cache"))
    def SetSwshCacheDirectory(self, value):
        self.swsh_cache_dir = value
        self.Modified()

    def _get_timesteps(self):
        logger.debug("Getting time range from data...")
        waveform_data = self._get_waveform_data()
        if len(waveform_data.RowData.keys()) == 0:
            return None
        ts = waveform_data.RowData["Time"]
        # XXX: why do we need linspace?
        return np.linspace(ts[0], ts[-1], len(ts))

    @smproperty.doublevector(
        name="TimestepValues",
        information_only="1",
        si_class="vtkSITimeStepsProperty",
    )
    def GetTimestepValues(self):
        timesteps = self._get_timesteps()
        return timesteps.tolist() if timesteps is not None else None

    def RequestInformation(self, request, inInfo, outInfo):
        logger.debug("Requesting information...")
        waveform_data_info = inInfo[0].GetInformationObject(0)

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

        return 1

    def RequestData(self, request, inInfo, outInfo):
        logger.info("Requesting data...")
        waveform_data = self._get_waveform_data()
        # grid_data = self._get_grid_data()
        output = dsa.WrapDataObject(vtkUniformGrid.GetData(outInfo))

        t = get_timestep(self, logger=logger)
        N = self.num_points_per_dim
        D = self.size

        dx = 2.0 * D / N
        N_y = N
        N_z = N
        output.SetDimensions(N, N_y, N_z)
        output.SetOrigin(-D, -D, -D)
        output.SetSpacing(dx, dx, dx)

        grid_params = SWSHParameters()
        grid_params.size = D
        grid_params.num_points = N
        grid_params.spin_weight = self.spin_weight
        grid_params.ell_max = self.ell_max

        swsh_grid = load_or_create_swsh_grid(grid_params, cache_dir=self.swsh_cache_dir)
        if swsh_grid is None:
            raise Exception('SWSH grid is None')

        adjust_params = GridAdjustParameters()
        adjust_params.scale_with_distance = self.scale_with_distance

        swsh_grid, spherical_grid = adjust_swsh_grid(swsh_grid, grid_params, adjust_params)

        # Compute scaled waveform phase on the grid
        phase = (t + self.time_shift) - spherical_grid.r

        # Compute quantity in the volume from the input waveform data
        waveform_timesteps = waveform_data.RowData["Time"]
        quantity = np.zeros(len(spherical_grid.r), dtype=np.complex)

        for i in range(self.modes_selection.GetNumberOfArrays()):
            mode_name = self.modes_selection.GetArrayName(i)
        for l in range(abs(grid_params.spin_weight), grid_params.ell_max + 1):
            for abs_m in range(0, l + 1):
                mode_name = get_mode_name(l, abs_m)
                quantity_mode = np.zeros(len(spherical_grid.r), dtype=np.complex)
                if not self.modes_selection.ArrayIsEnabled(mode_name):
                    continue
                for sign_m in (-1, 1):
                    m = abs_m * sign_m
                    dataset_name = "Y_l{}_m{}".format(l, m)
                    mode_profile = swsh_grid[:, LM_index(l, m, 0)]
                    waveform_mode_data = waveform_data.RowData[dataset_name]
                    if isinstance(waveform_mode_data, dsa.VTKNoneArray):
                        logger.warning(
                            f"Dataset '{dataset_name}' for mode {(l, m)} not"
                            " available in waveform data, skipping."
                        )
                        continue

                    waveform_mode_data = (
                        waveform_mode_data[:, 0] + 1j * waveform_mode_data[:, 1]
                    )
                    mode_data = np.interp(
                        phase,
                        waveform_timesteps,
                        waveform_mode_data,
                        left=0.0,
                        right=0.0,
                    )
                    quantity_mode += mode_data * mode_profile
                quantity += quantity_mode


        real_energy_flux = np.real(quantity)

        # Spherical clip
        real_energy_flux *= np.where(spherical_grid.r < grid_params.size*0.96, 1, 0)

        # Clip from below
        np.maximum(real_energy_flux, self.value_threshold, out=real_energy_flux)

        # Take log here instead of in ParaView
        np.log10(real_energy_flux, out=real_energy_flux)

        # Add entire sum to the output
        if self.component_selection.ArrayIsEnabled("Energy flux"):
            quantity_real_vtk = vtknp.numpy_to_vtk(real_energy_flux, deep=True)
            quantity_real_vtk.SetName("Energy flux")
            output.GetPointData().AddArray(quantity_real_vtk)

        return 1
