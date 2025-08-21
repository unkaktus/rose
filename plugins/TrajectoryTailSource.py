# Trajectory data Paraview reader

import logging

import sxs

import sys

import h5py
import numpy as np
from paraview import vtk
from paraview.util.vtkAlgorithm import smdomain, smhint, smproperty, smproxy
from paraview.vtk.util import numpy_support as vtknp
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.vtkCommonDataModel import vtkPolyData

logger = logging.getLogger(__name__)

def set_output_array(output, name, array):
    quantity_vtk = vtknp.numpy_to_vtk(array, deep=True)
    quantity_vtk.SetName(name)
    output.GetPointData().AddArray(quantity_vtk)

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

@smproxy.source(name="TrajectoryTailSource", label="Trajectory Tail Source")
class TrajectoryTailSource(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=0, nOutputPorts=1, outputType="vtkPolyData"
        )
        self.simulation_name = ""
        self.time = None
        self.coors = None

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

    @smproperty.stringvector(name="SimulationName", default_values="")
    def SetSimulationName(self, value):
        self.simulation_name = value
        self.Modified()

    @smproperty.doublevector(name="TotalMass", default_values=1.0)
    def SetTotalMass(self, value):
        self.total_mass = value
        self.Modified()

    @smproperty.stringvector(name="HorizonName", default_values="A")
    def SetSubfile(self, value):
        self.horizon_name = value
        self.Modified()

    @smproperty.doublevector(name="MinimumAge", default_values=0.0)
    def SetMinimumAge(self, value):
        self.minimum_age = value
        self.Modified()

    @smproperty.doublevector(name="MaximumAge", default_values=800.0)
    def SetMaximumAge(self, value):
        self.maximum_age = value
        self.Modified()

    def RequestInformation(self, request, inInfo, outInfo):
        logger.debug("Requesting information...")
        info = outInfo.GetInformationObject(0)

        if self.simulation_name == "":
            return -1

        sim = sxs.load(self.simulation_name, download=False)
        horizons = sim.horizons
        
        horizon = None
        if self.horizon_name == "A":
            horizon = horizons.A
        if self.horizon_name == "B":
            horizon = horizons.B
        if self.horizon_name == "C":
            horizon = horizons.C

        self.time = horizon.time * self.total_mass

        self.coords = np.column_stack((
            horizon.coord_center_inertial[:,0],
            horizon.coord_center_inertial[:,1],
            horizon.coord_center_inertial[:,2],
        ))
        self.coords *= self.total_mass

        set_timesteps(self, self._get_timesteps())
        return 1


    def RequestData(self, request, inInfo, outInfo):
        output = dsa.WrapDataObject(vtkPolyData.GetData(outInfo))

        if self.simulation_name == "":
            return -1

        current_time = get_timestep(self)
        age = current_time - self.time
        # Filter the point data
        mask = np.logical_and(self.minimum_age < age, age < self.maximum_age)

        # Construct a line of points
        points_vtk = vtk.vtkPoints()
        # Each ID is composed of (1) the order of the point in the line and (2)
        # the index in the `vtkPoints` constructed above
        line_vtk = vtk.vtkPolyLine()
        point_ids = line_vtk.GetPointIds()
        point_ids.SetNumberOfIds(len(self.coords[mask]))


        for i, point in enumerate(self.coords[mask]):
            points_vtk.InsertPoint(i, *point)
            point_ids.SetId(i, i)
        output.SetPoints(points_vtk)
        # Set the line ordering as "cell data"
        output.Allocate(1, 1)
        output.InsertNextCell(line_vtk.GetCellType(), line_vtk.GetPointIds())

        # Add data to the output
        set_output_array(output, name="Time", array=self.time[mask])
        set_output_array(output, name="Age", array=age[mask])

        return 1
