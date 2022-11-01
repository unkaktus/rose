# WaveformToVolume Paraview filter

import spherical
import quaternionic
import numpy as np
import os
import hashlib
import time
import logging
import sys
import tempfile

import rich.progress

from vtkmodules.vtkCommonDataModel import vtkUniformGrid
from vtkmodules.vtkCommonCore import vtkDataArraySelection
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.numpy_interface import dataset_adapter as dsa
from paraview.vtk.util import numpy_support as vtknp
from paraview.util.vtkAlgorithm import smdomain, smproperty, smproxy
from paraview import util


logger = logging.getLogger(__name__)
tmp_dir = tempfile.mkdtemp(prefix='rose-')


def new_progress():
    progress = rich.progress.Progress(
        rich.progress.TextColumn(
            "[progress.description]{task.description}"
        ),
        rich.progress.SpinnerColumn(
            spinner_name="simpleDots", finished_text="... done."
        ),
        rich.progress.TimeElapsedColumn(),
    )
    if progress is None:
        logger.fatal("rich progress is None")
        return
    return progress


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
    key = f'{p.size}|{p.num_points}|{p.spin_weight}|{p.ell_max}|{p.clip_y_normal}|{p.clip_z_normal}'.encode(
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
    progress = new_progress()

    task_id = progress.add_task("Computing SWSH grid", total=1)

    cartesian_grid = CartesianGrid(p)
    spherical_grid = cartesian_grid.spherical()

    angles = quaternionic.array.from_spherical_coordinates(
        spherical_grid.theta, spherical_grid.phi)
    swsh_grid = spherical.Wigner(p.ell_max).sYlm(s=p.spin_weight, R=angles)
    progress.update(task_id, completed=1)
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
    radial_scale: float
    one_over_r_scaling: bool


def adjust_swsh_grid(swsh_grid, grid_params: SWSHParameters, params: GridAdjustParameters):
    r = CartesianGrid(grid_params).spherical().r
    screen = activation(
        x=r - params.activation_offset,
        width=params.activation_width,
    ) * deactivation(x=r,
                     width=params.deactivation_width,
                     outer=grid_params.size)
    swsh_grid *= screen.reshape(screen.shape + (1,))
    # Apply radial scale
    r *= params.radial_scale
    if params.one_over_r_scaling:
        swsh_grid /= (r + 1.0e-30).reshape(r.shape + (1,))
    return swsh_grid


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


# Reproduces `spherical_functions.LM_index` so we don't need to import the
# `spherical_functions` module when using a cached SWSH grid
def LM_index(ell, m, ell_min):
    return ell * (ell + 1) - ell_min**2 + m


@smproxy.filter(label="Waveform To Volume")
@smproperty.input(name="WaveformData", port_index=0)
@smdomain.datatype(dataTypes=["vtkTable"])
class WaveformToVolume(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self,
            nInputPorts=1,
            nOutputPorts=1,
            # Choosing `vtkUniformGrid` for the output for the following reasons:
            # - `vtkRectilinearGrid` doesn't support volume rendering
            #   (in Paraview v5.7.0 at least)
            # - The unstructured grids don't support the 'GPU Based'
            #   volume rendering mode, which can do shading and looks nice
            outputType="vtkUniformGrid",
        )
        self.modes_selection = vtkDataArraySelection()
        # TODO: We should really retrieve the available modes from the input
        # info in `RequestInformation`, but the `WAVEFORM_MODES` keys is not
        # propagating downstream for some reason...
        # modes_arrays_key = WaveformDataReader.MODES_ARRAYS_KEY
        # if waveform_data_info.Has(modes_arrays_key):
        #     for i in range(waveform_data_info.Length(modes_arrays_key)):
        #         self.modes_selection.AddArray(waveform_data_info.Get(
        #             modes_arrays_key, i))
        for l in range(2, 5 + 1):
            for m in range(0, l + 1):
                self.modes_selection.AddArray(get_mode_name(l, m))
        self.modes_selection.AddObserver(
            "ModifiedEvent", create_modified_callback(self)
        )
        self.polarizations_selection = vtkDataArraySelection()
        self.polarizations_selection.AddArray("Plus")
        self.polarizations_selection.AddArray("Cross")
        self.polarizations_selection.AddObserver(
            "ModifiedEvent", create_modified_callback(self)
        )

    def FillInputPortInformation(self, port, info):
        # When using multiple inputs we (may) have to set their data types here
        # info.Set(self.INPUT_REQUIRED_DATA_TYPE(),
        #          'vtkTable' if port == 0 else 'vtkUniformGrid')
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

    @smproperty.intvector(name="NormalizeEachMode", default_values=False)
    def SetNormalizeEachMode(self, value):
        self.normalize_each_mode = value
        self.Modified()

    @smproperty.dataarrayselection(name="Polarizations")
    def GetPolarizations(self):
        return self.polarizations_selection

    # Not needed when using SwshGrid input
    @smproperty.doublevector(name="Size", default_values=100)
    def SetSize(self, value):
        self.size = value
        self.Modified()

    @smproperty.intvector(name="SpatialResolution", default_values=100)
    def SetSpatialResolution(self, value):
        self.num_points_per_dim = value
        self.Modified()

    @smproperty.intvector(name="KeepEveryNthTimestep", default_values=1)
    def SetKeepEveryNthTimestep(self, value):
        self.keep_every_n_timestep = value
        self.Modified()

    @smproperty.intvector(name="EllMax", default_values=2)
    def SetEllMax(self, value):
        self.ell_max = value
        self.Modified()

    @smproperty.intvector(name="SpinWeight", default_values=-2)
    def SetSpinWeight(self, value):
        self.spin_weight = value
        self.Modified()

    @smproperty.doublevector(name="RadialScale", default_values=10)
    def SetRadialScale(self, value):
        self.radial_scale = value
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

    @smproperty.intvector(name="InvertRotationDirection", default_values=False)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetInvertRotationDirection(self, value):
        self.invert_rotation_direction = value
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

    @smproperty.stringvector(name="SwshCacheDirectory", default_values=os.path.join(tmp_dir, "swsh_cache"))
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

        # For the `vtkUniformGrid` output we need to provide extents
        # so that it gets rendered at all.
        # When using the SwshGrid input we can retrieve them from the
        # information object and pass them on.
        # grid_extents = grid_info.Get(self.GetExecutive().WHOLE_EXTENT())
        N = self.num_points_per_dim
        N_y = N // 2 if self.clip_y_normal else N
        N_z = N // 2 if self.clip_z_normal else N
        grid_extents = [0, N - 1, 0, N_y - 1, 0, N_z - 1]
        util.SetOutputWholeExtent(self, grid_extents)

        # This needs the time data from the waveform file, so we may have to
        # set the `TIME_RANGE` and `TIME_STEPS` already in the
        # WaveformDataReader.
        set_timesteps(self, self._get_timesteps(), logger=logger)

        # logger.debug("Information object: {}".format(info))
        return 1

    def RequestData(self, request, inInfo, outInfo):
        logger.info("Requesting data...")
        waveform_data = self._get_waveform_data()
        # grid_data = self._get_grid_data()
        output = dsa.WrapDataObject(vtkUniformGrid.GetData(outInfo))

        t = get_timestep(self, logger=logger)
        N = self.num_points_per_dim
        D = self.size

        # We may have to forward the grid data here when using SwshGrid input
        # output.SetDimensions(*grid_data.GetDimensions())
        # output.SetOrigin(*grid_data.GetOrigin())
        # output.SetSpacing(*grid_data.GetSpacing())
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

        # Apply activation, radial scale etc.
        adjust_params = GridAdjustParameters()

        adjust_params.radial_scale = self.radial_scale
        adjust_params.activation_offset = self.activation_offset
        adjust_params.activation_width = self.activation_width
        adjust_params.deactivation_width = self.deactivation_width
        adjust_params.one_over_r_scaling = self.one_over_r_scaling

        swsh_grid = adjust_swsh_grid(swsh_grid, grid_params, adjust_params)

        spherical_grid = CartesianGrid(grid_params).spherical()

        logger.info(f"Computing volume data at t={t}...")
        start_time = time.time()

        # Compute scaled waveform phase on the grid
        # r = vtknp.vtk_to_numpy(grid_data.GetPointData()['RadialCoordinate'])
        phase = t - spherical_grid.r + self.activation_offset * self.radial_scale

        # Invert rotation direction
        rotation_direction = -1.0 if self.invert_rotation_direction else 1.0

        # Compute strain in the volume from the input waveform data
        skip_timesteps = self.keep_every_n_timestep
        waveform_timesteps = waveform_data.RowData["Time"][::skip_timesteps]
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
                    dataset_name = "Y_l{}_m{}".format(l, m)
                    mode_profile = swsh_grid[:, LM_index(l, m, 0)]
                    # mode_profile = vtknp.vtk_to_numpy(grid_data.GetPointData()[dataset_name])
                    waveform_mode_data = waveform_data.RowData[dataset_name][
                        ::skip_timesteps
                    ]
                    if isinstance(waveform_mode_data, dsa.VTKNoneArray):
                        logger.warning(
                            f"Dataset '{dataset_name}' for mode {(l, m)} not"
                            " available in waveform data, skipping."
                        )
                        continue
                    # TODO: Make sure inverting the rotation direction like this
                    # is correct.
                    waveform_mode_data = (
                        waveform_mode_data[:, 0]
                        + rotation_direction * 1j * waveform_mode_data[:, 1]
                    )
                    if self.normalize_each_mode:
                        waveform_mode_data /= np.max(
                            np.abs(waveform_mode_data))
                    mode_data = np.interp(
                        phase,
                        waveform_timesteps,
                        waveform_mode_data,
                        left=0.0,
                        right=0.0,
                    )
                    strain_mode += mode_data * mode_profile
                strain += strain_mode
                # Expose individual modes in output
                if self.store_individual_modes:
                    if self.polarizations_selection.ArrayIsEnabled("Plus"):
                        strain_mode_real_vtk = vtknp.numpy_to_vtk(
                            np.real(strain_mode), deep=True
                        )
                        strain_mode_real_vtk.SetName(mode_name + " Plus")
                        output.GetPointData().AddArray(strain_mode_real_vtk)
                    if self.polarizations_selection.ArrayIsEnabled("Cross"):
                        strain_mode_imag_vtk = vtknp.numpy_to_vtk(
                            np.imag(strain_mode), deep=True
                        )
                        strain_mode_imag_vtk.SetName(mode_name + " Cross")
                        output.GetPointData().AddArray(strain_mode_imag_vtk)

        if self.polarizations_selection.ArrayIsEnabled("Plus"):
            strain_real_vtk = vtknp.numpy_to_vtk(np.real(strain), deep=True)
            strain_real_vtk.SetName("Plus strain")
            output.GetPointData().AddArray(strain_real_vtk)

        if self.polarizations_selection.ArrayIsEnabled("Cross"):
            strain_imag_vtk = vtknp.numpy_to_vtk(np.imag(strain), deep=True)
            strain_imag_vtk.SetName("Cross strain")
            output.GetPointData().AddArray(strain_imag_vtk)

        logger.info(
            f"Volume data computed in {time.time() - start_time:.3f}s.")
        return 1
