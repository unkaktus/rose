## rose

`rose` is a framework to visualize gravitational-wave radiation using ParaView.
It is based on plugins from [`gwpv`](https://github.com/nilsvu/gwpv) written by [Nils Vu](https://github.com/nilsvu).

The name `rose` owes to the rose-like shapes of BBH gravitaional emission using `inferno` colormap.

![The rose](rose.png)

### Architecture
- ParaView plugins to load waveform and horizon data
- ParaView server packaged into a container together with the plugins and their dependencies
- Python script `rose-state.py` to render ParaView state into frames
- TODO: Slurm SBATCH job script to submit `rose-state.py` for rendering different frame ranges in parallel

### Pipeline

### Installation
1. Load (or install) Apptainer (previously known as Singularity).
For example:

```shell
$ source /home/SPACK2023/share/spack/setup-env.sh
$ module avail apptainer
$ module load apptainer-1.0.3-gcc-12.2.0-aojy6ca
```

2. Create a directory to store executable container files:
```shell
$ mkdir $HOME/apptainers
$ export PATH=$HOME/apptainers:$PATH
```
You might want to add it to `PATH` into your `.bashrc` file permanently.

3. There, download latest `rose` container file from GitHub:
```shell
$ cd ~/apptainers
$ wget https://github.com/unkaktus/rose/releases/download/v0.0.1/rose.sif
$ chmod +x rose.sif
```

### Usage

1. Run the container on some node:
```shell
$ srun -N1 -n1 --exclusive --pty rose.sif
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

6. Enjoy loading files from the cluster using `WaveformDataReader`, `TrajectoryDataReader`, and then applying filters `WaveformToVolume` or `TrajectoryTail`!
