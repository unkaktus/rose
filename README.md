## rose

![The rose](rose.png)

`rose` is a framework to visualize gravitational-wave radiation using ParaView.

It is based on plugins from [`gwpv`](https://github.com/nilsvu/gwpv) written by [Nils Vu](https://github.com/nilsvu).

The name `rose` comes from the flower shapes of gravitaional emission of compact binary coalescences.

> Detrás de los zarzales salvajes de tu pecho \
> Hay una rosa que deslumbrará todo el jardín

_Ruido_ - La Prohibida

### Workflow
1. Start ParaView server from the container
2. Connect to the server in ParaView desktop client
3. Setup the visualization using `rose` plugins
4. Save the state to a file
5. Copy the file to the cluster to render
6. Run the rendering job to render frames for the state
7. Run the overplotting job to create frames with the data, text, legends and logos
8. Run the job to combine the frames to a video file


### Installation

1. Install MambaForge locally and on the cluster.
Run the code below and follow the instructions, activating `base` environment:

```shell
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh
```

_In case you don't have internet access on the remote:_ use `mitten` instead of `ssh` to login there.
Install `mitten` locally via:
```shell
mamba install -c https://mamba.unkaktus.art/ mitten
```

2. Locally, install ParaView, spanner:
```shell
mamba install -c https://mamba.unkaktus.art/ paraview=5.11 spanner
```

3. On the cluster, install apptainer and spanner:
```shell
mamba install -c https://mamba.unkaktus.art/ apptainer spanner
```

4. Create a directory for containers and download the latest `rose` container file from GitHub:
```shell
mkdir -p $HOME/apptainers
curl -L -o $HOME/apptainers/rose-v2.0.0.sif https://github.com/unkaktus/rose/releases/download/v2.0.0/rose-v2.0.0.sif
```


### Running on an HPC cluster in interactive mode

1. Run the ParaView server from the container, adjusting the directories you want to mount inside the container appropriately:
```shell
srun -J rose --pty apptainer exec --bind /scratch:/scratch --bind /work:/work ~/apptainers/rose-v2.0.0.sif pvserver
```

2. Locally, start ParaView by running `paraview` from the terminal.

3. Setup port forwarding to the remote ParaView server job, adjusting the SSH machine name and name of the job:
```shell
spanner port-forward -p 11111 -m machine-name rose
```

4. In the ParaView window, connect to remote ParaView server - `Connect` button.
Then, specify `localhost:11111` as the address of the server.

6. Enjoy loading files from the cluster using `EnergyFluxVolumeReader`, `StrainVolumeReader`, and `TrajectoryTailReader`!



### Running parallel state rendering on an HPC cluster

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
