import logging

logger = logging.getLogger(__name__)

import numpy as np
import h5py
from paraview.util.vtkAlgorithm import smdomain, smhint, smproperty, smproxy
from paraview.vtk.util import keys as vtkkeys
from paraview.vtk.util import numpy_support as vtknp
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.vtkCommonDataModel import vtkTable

import scri
from spherical_functions import LM_index


@smproxy.reader(
    name="EnergyFluxReader",
    label="Energy Flux Reader",
    extensions="h5",
    file_description="HDF5 files",
)
class EnergyFluxReader(VTKPythonAlgorithmBase):
    WAVEFORM_MODES_KEY = vtkkeys.MakeKey(
        vtkkeys.StringVectorKey, "WAVEFORM_MODES", "EnergyFluxReader"
    )

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=0, nOutputPorts=1, outputType="vtkTable"
        )
        self._filename = None
        self._subfile = None
        self.mode_names = []

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

    def RequestInformation(self, request, inInfo, outInfo):
        info = outInfo.GetInformationObject(0)
        
        return 1

    def RequestData(self, request, inInfo, outInfo):
        output = dsa.WrapDataObject(vtkTable.GetData(outInfo))

        if (
            self._filename is not None
            and self._subfile is not None
        ):
            abd = scri.SpEC.create_abd_from_h5("SXS", h = f"{self._filename}/{self._subfile}")
            
            # Get the timesteps
            t = abd.t

            # Add timesteps column
            column_time = vtknp.numpy_to_vtk(t, deep=False)
            column_time.SetName("Time")
            output.AddColumn(column_time)

            ## XXX set it manually for now
            self.spin_weight = 0
            self.ell_max = 8


            # Cache calculated energy flux
            energy_flux = np.array([])
            energy_flux_filename = f'{self._filename}.energy_flux.npy'

            try:
                energy_flux = np.load(energy_flux_filename)
            except:
                # Calculate the energy flux modes coefficients
                energy_flux = np.array(abd.sigma.dot * abd.sigma.dot.bar)
                np.save(energy_flux_filename, energy_flux)

            
            # Add modes
            for l in range(np.abs(self.spin_weight), self.ell_max + 1):
                for m in range(0, l + 1):
                    abs_m = abs(m)
                    for sign_m in (-1, 1):
                        m = sign_m * abs_m
                        mode_name = f'Y_l{l}_m{m}'

                        # Extract coefficients for the mode
                        mode_energy_flux = np.array(energy_flux[:, LM_index(l,m,0)])

                        # Add energy flux column
                        mode_energy_flux_column = np.column_stack((mode_energy_flux.real,mode_energy_flux.imag))
                        column_mode = vtknp.numpy_to_vtk(mode_energy_flux_column, deep=False)
                        column_mode.SetName(mode_name)
                        output.AddColumn(column_mode)

        return 1
