# Trajectory tail Paraview filter

import logging

import numpy as np
from paraview.util.vtkAlgorithm import smdomain, smproperty, smproxy
from paraview.vtk.util import numpy_support as vtknp
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.vtkCommonDataModel import vtkPolyData

logger = logging.getLogger(__name__)


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

@smproxy.filter(label="Trajectory Tail")
@smproperty.input(name="TrajectoryData", port_index=0)
@smdomain.datatype(dataTypes=["vtkPolyData"])
class TrajectoryTail(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=1, nOutputPorts=1, outputType="vtkPolyData"
        )

    def FillInputPortInformation(self, port, info):
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData")

    def _get_trajectory_data(self):
        return dsa.WrapDataObject(self.GetInputDataObject(0, 0))

    def _get_timesteps(self):
        logger.debug("Getting time range from data...")
        trajectory_data = self._get_trajectory_data()
        point_times = trajectory_data.PointData["Time"]
        # Using a few timesteps within the data range so we can animate through
        # them in the GUI
        return np.linspace(point_times[0], point_times[-1], 100)

    @smproperty.doublevector(
        name="TimestepValues",
        information_only="1",
        si_class="vtkSITimeStepsProperty",
    )
    def GetTimestepValues(self):
        return self._get_timesteps().tolist()

    def RequestInformation(self, request, inInfo, outInfo):
        logger.debug("Requesting information...")
        # This needs the time data from the trajectory file, so we may have to
        # set the `TIME_RANGE` and `TIME_STEPS` already in the
        # TrajectoryDataReader.
        set_timesteps(self, self._get_timesteps(), logger=logger)
        return 1

    def RequestData(self, request, inInfo, outInfo):
        logger.debug("Tail: Requesting data...")
        logger.debug(f"{inInfo}")

        input = self.GetInputDataObject(0, 0)
        trajectory_data = dsa.WrapDataObject(input)

        output = dsa.WrapDataObject(vtkPolyData.GetData(outInfo))

        # Shallow-copy input trajectory data to output
        output.ShallowCopy(input)

        # Retrieve current time
        time = get_timestep(self, logger=logger)

        # Add age data to the points
        age = time - trajectory_data.PointData["Time"]
        age_vtk = vtknp.numpy_to_vtk(age, deep=True)
        age_vtk.SetName("Age")
        output.GetPointData().AddArray(age_vtk)

        logger.debug("Tail: Requesting data done")

        return 1
