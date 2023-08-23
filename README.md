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

### Requirements
* ParaView 5.11.0 on your desktop
* [`spanner`](https://github.com/unkaktus/spanner) for running render in parallel


### Installation on an HPC cluster
1. Load or install Apptainer. It is likely it is modules on your system.

For example:

```shell
source /home/SPACK2023/share/spack/setup-env.sh
module avail apptainer
module load apptainer-1.0.3-gcc-12.2.0-aojy6ca
```

2. Create a directory to store executable container files:
```shell
mkdir $HOME/apptainers
```

3. There, download latest `rose` container file from GitHub:
```shell
cd ~/apptainers
wget https://github.com/unkaktus/rose/releases/download/v2.0.0/rose-v2.0.0.sif
```

### Running on an HPC cluster in interactive mode

1. Run the container on some node:
```shell
srun -J rose-interactive -N1 -n1 --exclusive --pty apptainer shell --bind /scratch:/scratch --bind /work:/work ~/apptainers/rose-v2.0.0.sif
```
2. Note the hostname of the node you are running at:
```shell
hostname
```
3. Run ParaView server inside the container:
```shell
pvserver
```

4. Create a new SSH connection with local port forwarding of the port 11111 to the `pvserver` on the node it is running on. Here, replace `node-hostname` with the node hostname from Step 2, as well as your cluster hostname and username.
```shell
ssh -L 11111:node-hostname:11111 username@cluster
```
5. In local ParaView, connect to remote ParaView server - `Connect` button.
Then, specify `localhost:11111` as the address of the server.

6. Enjoy loading files from the cluster using `EnergyFluxVolumeReader`, `StrainVolumeReader`, and `TrajectoryTailReader`!

### Running parallel state rendering on an HPC cluster

0. Make sure you installed `spanner` (see above).

1. Copy `render/begin.toml` into the directory with your state file (`.pvsm`).

2. Adjust the number of nodes, path to `rose.sif` container, and other parameters in `begin.toml`.

3. Begin the rendering job on your state using `spanner`:
```shell
spanner begin state.pvsm
```
4. You can track the process using `spanner`, e.g.:
```shell
spanner logs -f rose_state
```

5. The job will create a directory with the same name as the state file, and output the frames there.

6. You can create a video from these frames by running `make_movie.sh` in the frame directory:
```shell
apptainer exec --bind /scratch:/scratch --bind /work:/work ~/apptainers/rose-v2.0.0.sif make_movie.sh
```

This will create a video file with the name of the current directory.

### Running parallel frame postprocessing on an HPC cluster

0. Make sure you installed `spanner` (see above).

1. Copy `render/begin-overplot.toml` into the directory with your state file (`.pvsm`).

2. Adjust the number of nodes, path to `rose.sif` container, and other parameters in `begin-overplot.toml`.

3. Begin the rendering job on your state using `spanner`:
```shell
spanner begin -f begin-overplot.toml state.pvsm
```
4. You can track the process using `spanner`, e.g.:
```shell
spanner logs -f rose-overplot_state
```

5. The job will plot legends and a colorbar on top of the frames in the frame directory, and save
postprocessed frames with `frame-overplotted` prefix.

6. You can create a video from these frames by running `make_movie.sh` in the frame directory using the correct frame name prefix:
```shell
apptainer exec --bind /scratch:/scratch --bind /work:/work ~/apptainers/rose-v2.0.0.sif make_movie.sh frame-overplotted
```

It will create a video file with the name of the current directory with the frame prefix included.

This is your final movie, enjoy! 💫
