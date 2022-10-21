FROM python:3.9.5

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -yqq update \
 && apt-get -yqq upgrade \
 && apt-get -yqq install --no-install-recommends \
      curl \
      ffmpeg \
      git \
      vim-tiny \
      make \
      wget \
      mpich \
      libhdf5-dev \
  && apt-get clean

# Install ParaView
RUN wget -O paraview.tar.gz "https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.10&type=binary&os=Linux&downloadFile=ParaView-5.10.1-osmesa-MPI-Linux-Python3.9-x86_64.tar.gz" \
  && tar -xzf paraview.tar.gz \
  && rm paraview.tar.gz \
  && mv ParaView-* /opt/paraview

ENV PYTHONPATH="/opt/paraview/lib/python3.9/site-packages:/usr/local/lib/python3.9/dist-packages:/usr/local/lib/python3.9/site-packages" \
  PATH="/opt/paraview/bin:$PATH"
ENV LC_ALL="C.UTF-8"

# Install Python dependencies
RUN --mount=type=cache,target=/root/.cache pip install \
    astropy \
    h5py>=3.0.0 \
    numpy \
    rich \
    scipy \
    spherical \
    quaternionic \
    ;

WORKDIR /opt/rose
COPY . .
ENV PV_PLUGIN_PATH=/opt/rose/plugins
ENTRYPOINT ["/bin/bash"]