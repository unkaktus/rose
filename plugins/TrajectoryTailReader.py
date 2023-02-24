# Trajectory data Paraview reader

import logging

import h5py
import numpy as np
from paraview import vtk
from paraview.util.vtkAlgorithm import smdomain, smhint, smproperty, smproxy
from paraview.vtk.util import numpy_support as vtknp
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.vtkCommonDataModel import vtkPolyData

logger = logging.getLogger(__name__)


def get_timestep(algorithm):
    logger.debug("Getting current timestep...")
    executive = algorithm.GetExecutive()
    outInfo = executive.GetOutputInformation(0)
    if not outInfo.Has(executive.UPDATE_TIME_STEP()):
        logger.debug("No `UPDATE_TIME_STEP` found. Return 0.")
        return 0.0
    timestep = outInfo.Get(executive.UPDATE_TIME_STEP())
    logger.debug(f"Found `UPDATE_TIME_STEP`: {timestep}")
    return timestep

def set_timesteps(algorithm, timesteps):
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


@smproxy.reader(
    label="Trajectory Tail Reader",
    extensions="h5",
    file_description="HDF5 files",
)
class TrajectoryTailReader(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=0, nOutputPorts=1, outputType="vtkPolyData"
        )

    def _get_timesteps(self):
        return np.linspace(self.time[0], self.time[-1], len(self.time))

    @smproperty.doublevector(
        name="TimestepValues",
        information_only="1",
        si_class="vtkSITimeStepsProperty",
    )
    def GetTimestepValues(self):
        timesteps = self._get_timesteps()
        return timesteps.tolist() if timesteps is not None else None


    @smproperty.stringvector(name="File")
    @smdomain.filelist()
    @smhint.filechooser(extensions="h5", file_description="HDF5 files")
    def SetFile(self, value):
        self._filename = value
        self.Modified()

    @smproperty.stringvector(name="Subfile", default_values="AhB.dir")
    def SetSubfile(self, value):
        self._subfile = value
        self.Modified()

    @smproperty.stringvector(
        name="CoordinatesDataset", default_values="CoordCenterInertial.dat"
    )
    def SetCoordinatesDataset(self, value):
        self._coords_dataset = value
        self.Modified()

    @smproperty.doublevector(name="RadialScale", default_values=1.0)
    def SetRadialScale(self, value):
        self._radial_scale = value
        self.Modified()


    def RequestInformation(self, request, inInfo, outInfo):
        logger.debug("Requesting information...")
        info = outInfo.GetInformationObject(0)

        trajectory_file = h5py.File(self._filename, "r")
        subfile = trajectory_file[self._subfile]
        coords = np.array(subfile[self._coords_dataset])

        coords[:, 1:] *= self._radial_scale

        self.time = coords[:,0]
        self.coords = coords[:,1:]
        # This needs the time data from the trajectory file, so we may have to
        # set the `TIME_RANGE` and `TIME_STEPS` already in the
        # TrajectoryDataReader.
        set_timesteps(self, self._get_timesteps())
        return 1


    def RequestData(self, request, inInfo, outInfo):
        output = dsa.WrapDataObject(vtkPolyData.GetData(outInfo))

        # Construct a line of points
        points_vtk = vtk.vtkPoints()
        # Each ID is composed of (1) the order of the point in the line and (2)
        # the index in the `vtkPoints` constructed above
        line_vtk = vtk.vtkPolyLine()
        point_ids = line_vtk.GetPointIds()
        point_ids.SetNumberOfIds(len(self.coords))
        for i, point in enumerate(self.coords):
            points_vtk.InsertPoint(i, *point)
            point_ids.SetId(i, i)
        output.SetPoints(points_vtk)
        # Set the line ordering as "cell data"
        output.Allocate(1, 1)
        output.InsertNextCell(line_vtk.GetCellType(), line_vtk.GetPointIds())

        # Add time data to the points
        time_vtk = vtknp.numpy_to_vtk(self.time)
        time_vtk.SetName("Time")
        output.GetPointData().AddArray(time_vtk)
        
        # Retrieve current time
        current_time = get_timestep(self)
     
        # Add age data to the points
        age = current_time - self.time
        age_vtk = vtknp.numpy_to_vtk(age, deep=True)
        age_vtk.SetName("Age")
        output.GetPointData().AddArray(age_vtk)

        return 1
