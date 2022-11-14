#!/bin/bash

conda create -n af python=3.7
conda activate af
conda install -y -c conda-forge openmm==7.5.1 cudnn==8.2.1.32 cudatoolkit==11.0.3 pdbfixer==1.7
conda install -y -c bioconda hmmer==3.3.2 hhsuite==3.3.0 kalign2==2.04
pip install absl-py==0.13.0 biopython==1.79 chex==0.0.7 dm-haiku==0.0.4 dm-tree==0.1.6 immutabledict==2.0.0 jax==0.2.14 ml-collections==0.1.0 numpy==1.19.5 scipy==1.7.0 tensorflow==2.5.0 pandas==1.3.4 tensorflow-cpu==2.5.0 future matplotlib python-igraph csb psutil paramiko scikit-learn scipy pickle5
pip install --upgrade jax==0.2.14 jaxlib==0.1.69+cuda111 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
pip install git+https://github.com/deepmind/alphafold.git@v2.1.2
path=$(python -c 'import site; print(site.getsitepackages()[0])')
cd $path
curl https://raw.githubusercontent.com/deepmind/alphafold/main/docker/openmm.patch -o /tmp/openmm.patch
patch -p0 < /tmp/openmm.patch
