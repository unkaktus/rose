import spherical
import quaternionic
import numpy as np
import os
import hashlib
import logging
import tempfile
import h5py


from vtkmodules.vtkCommonDataModel import vtkUniformGrid
from vtkmodules.vtkCommonCore import vtkDataArraySelection
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.numpy_interface import dataset_adapter as dsa
from paraview.vtk.util import numpy_support as vtknp
from paraview.util.vtkAlgorithm import smdomain, smhint, smproperty, smproxy
from paraview import util


logger = logging.getLogger(__name__)
rose_cache_dir = os.environ.get("ROSE_CACHE_DIR", "")
if rose_cache_dir is None:
    rose_cache_dir = tempfile.mkdtemp(prefix='rose-')

class SWSHParameters:
    size: float = 0.0
    num_points: float = 0.0
    spin_weight: int = -2
    ell_max: int = 0
    clip_y_normal: bool = False
    clip_z_normal: bool = False


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
        Y = np.linspace(-p.size, 0, p.num_points //
                        2) if p.clip_y_normal else X
        Z = np.linspace(-p.size, 0, p.num_points //
                        2) if p.clip_z_normal else X
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
    key = f'{p.size}|{p.num_points}|{p.spin_weight}|{p.ell_max}|{p.clip_y_normal}|{p.clip_z_normal}'.encode('utf-8')
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

def set_output_array(output, name, array):
    quantity_vtk = vtknp.numpy_to_vtk(array, deep=True)
    quantity_vtk.SetName(name)
    output.GetPointData().AddArray(quantity_vtk)


def smoothstep(x):
    return np.where(x < 0, 0, np.where(x <= 1, 3 * x**2 - 2 * x**3, 1))


def activation(x, width):
    return smoothstep(x / width)


def deactivation(x, width, outer):
    return smoothstep((outer - x) / width)


class GridAdjustParameters():
    activation_offset: float
    activation_width: float
    deactivation_width: float
    one_over_r_scaling: bool


def adjust_swsh_grid(swsh_grid, grid_params: SWSHParameters, params: GridAdjustParameters):
    spherical_grid = CartesianGrid(grid_params).spherical()
    screen = activation(
        x = spherical_grid.r - params.activation_offset,
        width=params.activation_width,
    ) * deactivation(x = spherical_grid.r,
                     width=params.deactivation_width,
                     outer=grid_params.size)
    swsh_grid *= screen.reshape(screen.shape + (1,))
    if params.one_over_r_scaling:
        swsh_grid /= (spherical_grid.r + 1.0e-30).reshape(spherical_grid.r.shape + (1,))
    return swsh_grid, spherical_grid


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

def get_timestep(algorithm):
    executive = algorithm.GetExecutive()
    outInfo = executive.GetOutputInformation(0)
    if not outInfo.Has(executive.UPDATE_TIME_STEP()):
        return 0.0
    timestep = outInfo.Get(executive.UPDATE_TIME_STEP())
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
    return f"({l}, {abs_m}) Mode"


def LM_index(ell, m, ell_min):
    return ell * (ell + 1) - ell_min**2 + m

@smproxy.reader(
    name="StrainVolumeReader",
    label="Strain Volume Reader",
    extensions="h5",
    file_description="HDF5 files",
)
class StrainVolumeReader(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self,
            nInputPorts=0,
            nOutputPorts=1,
            outputType="vtkUniformGrid",
        )
        self._filename = None
        self._subfile = None
        self.mode_names = []
        self.modes_selection = vtkDataArraySelection()
        self.strain_modes = {}
        self.spin_weight = -2
        self.ell_max = 2
        self.swsh_cache_dir = os.path.join(rose_cache_dir, "swsh_cache")

        self.update_mode_selection()

        self.polarizations_selection = vtkDataArraySelection()
        self.polarizations_selection.AddArray("Plus")
        self.polarizations_selection.AddArray("Cross")
        self.polarizations_selection.AddObserver(
            "ModifiedEvent", create_modified_callback(self)
        )

    def update_mode_selection(self):
        for l in range(abs(self.spin_weight), self.ell_max + 1):
            for m in range(0, l + 1):
                self.modes_selection.AddArray(get_mode_name(l, m))
        self.modes_selection.AddObserver(
            "ModifiedEvent", create_modified_callback(self)
        )

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

    @smproperty.dataarrayselection(name="Modes")
    def GetModes(self):
        return self.modes_selection

    @smproperty.dataarrayselection(name="Polarizations")
    def GetPolarizations(self):
        return self.polarizations_selection

    @smproperty.doublevector(name="Size", default_values=100)
    def SetSize(self, value):
        self.size = value
        self.Modified()

    @smproperty.intvector(name="SpatialResolution", default_values=100)
    def SetSpatialResolution(self, value):
        self.num_points_per_dim = value
        self.Modified()

    @smproperty.intvector(name="EllMax", default_values=2)
    def SetEllMax(self, value):
        self.ell_max = value
        self.update_mode_selection()
        self.Modified()

    @smproperty.intvector(name="ClipYNormal", default_values=False)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetClipYNormal(self, value):
        self.clip_y_normal = value
        self.Modified()

    @smproperty.intvector(name="ClipZNormal", default_values=False)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetClipZNormal(self, value):
        self.clip_z_normal = value
        self.Modified()

    @smproperty.intvector(name="OneOverRScaling", default_values=False)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetOneOverRScaling(self, value):
        self.one_over_r_scaling = value
        self.Modified()

    @smproperty.doublevector(name="ActivationOffset", default_values=10)
    def SetActivationOffset(self, value):
        self.activation_offset = value
        self.Modified()

    @smproperty.doublevector(name="ActivationWidth", default_values=10)
    def SetActivationWidth(self, value):
        self.activation_width = value
        self.Modified()

    @smproperty.doublevector(name="DeactivationWidth", default_values=10)
    def SetDeactivationWidth(self, value):
        self.deactivation_width = value
        self.Modified()

    def _get_timesteps(self):
        ts = np.linspace(self.time[0], self.time[-1], len(self.time))
        return ts

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
        info = outInfo.GetInformationObject(0)

        if self._filename is None or self._subfile is None:
            return -1

        with h5py.File(self._filename, "r") as f:
            self.mode_names = list(
                map(
                    lambda dataset_name: dataset_name.replace(".dat", ""),
                    filter(
                        lambda dataset_name: dataset_name.startswith("Y_"),
                        f[self._subfile].keys(),
                    ),
                )
            )
            if len(self.mode_names) == 0:
                logger.warning(
                    "No waveform mode datasets (prefixed 'Y_') found in file"
                    f" '{self._filename}:{self._subfile}'."
                )
                return -1

            strain = f[self._subfile]
            # Take the time from (2,2) mode as the global time
            self.time = strain["Y_l2_m2.dat"][:, 0]

            for mode_name in self.mode_names:
                mode_data_column = strain[mode_name + ".dat"]
                self.strain_modes[mode_name] = np.array(
                    mode_data_column[:, 1] + 1j * mode_data_column[:, 2]
                    )

        N = self.num_points_per_dim
        N_y = N // 2 if self.clip_y_normal else N
        N_z = N // 2 if self.clip_z_normal else N
        grid_extents = [0, N - 1, 0, N_y - 1, 0, N_z - 1]
        util.SetOutputWholeExtent(self, grid_extents)

        set_timesteps(self, self._get_timesteps(), logger=logger)

        return 1

    def RequestData(self, request, inInfo, outInfo):
        logger.info("Requesting data...")
        output = dsa.WrapDataObject(vtkUniformGrid.GetData(outInfo))

        t = get_timestep(self)
        N = self.num_points_per_dim
        D = self.size

        dx = 2.0 * D / N
        N_y = N // 2 if self.clip_y_normal else N
        N_z = N // 2 if self.clip_z_normal else N
        output.SetDimensions(N, N_y, N_z)
        output.SetOrigin(-D, -D, -D)
        output.SetSpacing(dx, dx, dx)

        grid_params = SWSHParameters()
        grid_params.size = D
        grid_params.num_points = N
        grid_params.spin_weight = self.spin_weight
        grid_params.ell_max = self.ell_max
        grid_params.clip_y_normal = self.clip_y_normal
        grid_params.clip_z_normal = self.clip_z_normal

        swsh_grid = load_or_create_swsh_grid(grid_params, cache_dir=self.swsh_cache_dir)
        if swsh_grid is None:
            raise Exception('SWSH grid is None')

        adjust_params = GridAdjustParameters()
        adjust_params.activation_offset = self.activation_offset
        adjust_params.activation_width = self.activation_width
        adjust_params.deactivation_width = self.deactivation_width
        adjust_params.one_over_r_scaling = self.one_over_r_scaling

        swsh_grid, spherical_grid = adjust_swsh_grid(swsh_grid, grid_params, adjust_params)

        # Compute scaled waveform phase on the grid
        phase = t - spherical_grid.r + self.activation_offset

        # Compute strain in the volume
        strain = np.zeros(len(spherical_grid.r), dtype=np.complex)

        for i in range(self.modes_selection.GetNumberOfArrays()):
            mode_name = self.modes_selection.GetArrayName(i)
        for l in range(abs(grid_params.spin_weight), grid_params.ell_max + 1):
            for abs_m in range(0, l + 1):
                mode_name = get_mode_name(l, abs_m)
                strain_mode = np.zeros(len(spherical_grid.r), dtype=np.complex)
                if not self.modes_selection.ArrayIsEnabled(mode_name):
                    continue
                for sign_m in (-1, 1):
                    m = abs_m * sign_m
                    mode_profile = swsh_grid[:, LM_index(l, m, 0)]
                    mode_data = self.strain_modes[f"Y_l{l}_m{m}"]

                    mode_data = np.interp(
                        phase,
                        self.time,
                        mode_data,
                        left=0.0,
                        right=0.0,
                    )
                    strain_mode += mode_data * mode_profile
                strain += strain_mode

        if self.polarizations_selection.ArrayIsEnabled("Plus"):
            set_output_array(output, "Plus strain", np.real(strain))

        if self.polarizations_selection.ArrayIsEnabled("Cross"):
            set_output_array(output, "Cross strain", np.imag(strain))

        return 1
