## rose

![The rose](rose.png)

`rose` is a framework to visualize gravitational-wave radiation using ParaView.

It is based on [`gwpv`](https://github.com/nilsvu/gwpv) written by [Nils Vu](https://github.com/nilsvu). Check it out!

The name `rose` comes from the flower shapes of gravitaional emission of compact binary coalescences.

> Detrás de los zarzales salvajes de tu pecho \
> Hay una rosa que deslumbrará todo el jardín

_Ruido_ - La Prohibida

### Workflow
1. Create and activate Mamba/Conda environment with the Python version matching the one in ParaView
2. Install the dependencies
3. Start ParaView using `rose-pv` script to automatically load the plugins
4. Open the waveform files and setup the scene, and export the frames
5. Run the overplotting script to add legends, text, logos, and other data
6. Combine the frames into a video
5. Enjoy the results!


### Installation

1. Install MambaForge using the command below:

```shell
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh
```

_In case you don't have internet access on the remote:_ use [`mitten`](https://github.com/unkaktus/mitten) instead of `ssh` to pass your internet connection to the remote.

2. Find the exact Python version your ParaView has. Go to ParaView->About ParaView and note down the "Python Library Version". For example, my ParaView 5.13.1 has Python 3.10.13.

3. Create and activate the Mamba environment for `rose` with the matching Python version and dependencies:

```shell
mamba create -y -n rose python=3.10.13 numpy scipy psutil astropy h5py spherical scri spherical_functions
mamba activate rose
```
4. From the `rose` directory root, start ParaView via `rose-pv` script by specifying path to your ParaView binary, 

```shell
./rose-pv /path/to/paraview
```

or the application in case of macOS:

```shell
./rose-pv /Applications/ParaView-5.13.1.app
```

That will start the ParaView and load all `rose` plugins.

5. Now you are ready to open your waveform files using the appropriate reader. For example, `rhOverM_Asymptotic_GeometricUnits_CoM.h5` using `EnergyFluxVolumeReader`.
Note that at the momement `rose` supports only extrapolated waveforms in SXS catalog format, cf. Appendix A.3.1 of Boyle:2019kee.

### Running parallel state rendering on an HPC cluster

Running locally on your machine might be too slow or not fitting into RAM for high resolution. For that reason, you might want to run it on a more powerful cluster.

The idea here is to create the scene in ParaView locally, save it to a state file, and then render it on a remote cluster. As the paths to the data files might be different locally and on the cluster, one has to readjust them by editing the state file. As I didn't write yet the script to automatically swap the path, one has has to do it manually.

The below instructions are currently outdated.

1. Copy `render/render.begin` and your state file (`.pvsm`) to the cluster.

2. Adjust the number of nodes, path to `rose.sif` container, and other parameters in `render.begin`.

3. Begin the rendering job on your state using `spanner`:
```shell
spanner begin -f /path/to/render.begin state.pvsm
```
4. Check the status of the job:
```shell
spanner list
```

6. Track the process:
```shell
spanner logs -f rose_job_name
```

7. The job will create a directory with the same name as the state file, and output the frames there.

8. Create a video named `output.mp4` from the frames in the run directory, adding dark background:
```shell
apptainer exec --bind /scratch:/scratch --bind /work:/work ~/apptainers/rose-v2.0.0.sif \
     create_video.py --frames-dir rendered --background-color 0x161616 --output output.mp4
```

### Running parallel frame postprocessing on an HPC cluster

1. Copy `render/overplot/overplot.begin` to the cluster.

2. Adjust the number of nodes, path to `rose.sif` container, and other parameters in `overplot.begin`.

3. Begin the rendering job on your state:
```shell
spanner begin -f /path/to/overplot.begin state.pvsm
```

4. The job will plot legends and a colorbar on top of the frames in the frame directory, and save
postprocessed frames into `overplotted` subdirectory.

5. Create a video from these frames:
```shell
apptainer exec --bind /scratch:/scratch --bind /work:/work ~/apptainers/rose-v2.0.0.sif \
     create_video.py --frames-dir overplotted --output output-overplotted.mp4
```

This is your final movie, enjoy! 💫

_Aquí y ahora puede comenzar tu viaje._
