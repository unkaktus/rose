## rose

`rose` is a framework to visualize gravitational-wave radiation using ParaView.

It is based on plugins from [`gwpv`](https://github.com/nilsvu/gwpv) written by [Nils Vu](https://github.com/nilsvu).

The name `rose` comes from the rose-like shapes of gravitaional emission of binary blackholes.

![The rose](rose.png)

> The closer I get to you \
> A feeling comes over me (Me, too) \
> Falling closer, sweet as the gravity \
> The closer I get to you

[_The Closer I Get to You_ - Roberta Flack & Donny Hathaway](https://air.unkaktus.art/WWW92xgv94c)

### Architecture
- ParaView plugins to load waveform and horizon data
- ParaView server packaged into a container together with the plugins and their dependencies
- Python script `render/render_state.py` to render ParaView state into frames
- Slurm job script `render/job.sh` to submit `render/render_state.py` for rendering frame ranges in parallel

### Requirements
* ParaView 5.11.0 on your desktop

### Running locally via Docker

0. Get Docker via installing Docker Desktop or via https://get.docker.com/.

1. Run the Docker container, mounting your home folder as `/home` and grids cache directory `~/.cache` inside the container.
```shell
$ docker run --rm -ti -v "$HOME:/home" -v "$HOME/.cache:/cache" -e ROSE_CACHE_DIR="/cache" -p 12321:11111 unkaktus/rose:v2.0.0
$ pvserver
```

2. Once `pvserver` is running, you can connect to it from ParaView using address `localhost:12321`.

3. Enjoy loading files from the cluster using `WaveformDataReader`, `TrajectoryDataReader`, and then applying filters `WaveformToVolume` or `TrajectoryTail`!


### Installation on an HPC cluster
1. Load or install Apptainer. It is likely it is modules on your system.

For example:

```shell
$ source /home/SPACK2023/share/spack/setup-env.sh
$ module avail apptainer
$ module load apptainer-1.0.3-gcc-12.2.0-aojy6ca
```

2. Create a directory to store executable container files:
```shell
$ mkdir $HOME/apptainers
```

3. There, download latest `rose` container file from GitHub:
```shell
$ cd ~/apptainers
$ wget https://github.com/unkaktus/rose/releases/download/v2.0.0/rose-v2.0.0.sif
```

### Running on an HPC cluster

1. Run the container on some node:
```shell
$ srun -N1 -n1 --exclusive --pty apptainer shell --bind /scratch:/scratch --bind /home:/home ~/apptainers/rose-v2.0.0.sif
```
2. Note the hostname of the node you are running at:
```shell
$ hostname
```
3. Run ParaView server inside the container:
```shell
Apptainer> pvserver
```

4. Create a new SSH connection with local port forwarding of the port 11111 to the `pvserver` on the node it is running on. Here, replace `node-hostname` with the node hostname from Step 2, as well as your cluster hostname and username.
```shell
$ ssh -L 11111:node-hostname:11111 username@cluster.aei.mpg.de
```
5. In local ParaView, connect to remote ParaView server - `Connect` button.
Then, specify `localhost:11111` as the address of the server.

6. Enjoy loading files from the cluster using the plugins `EnergyFluxToVolume`, `WaveformDataReader`, and `WaveformToVolume`!
