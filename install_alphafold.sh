#!/bin/bash -i

eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda create -y -n af2 python=3.8
conda activate af2
conda install -y -c conda-forge openmm==7.5.1 cudatoolkit==11.2 cudnn==8.2.1.32 pdbfixer==1.7
conda install -y -c bioconda hmmer==3.3.2 hhsuite==3.3.0 kalign2==2.04
if conda install -y -c conda-forge cudatoolkit-dev; then
    echo 'Cudatoolkit-dev installed'
else
    echo 'ERROR: Conda installation of cudatoolkit failed. Please, install cudatoolkit afterwards.'
fi
pip install absl-py==1.0.0 biopython==1.79 chex==0.0.7 dm-haiku==0.0.9 dm-tree==0.1.6 immutabledict==2.0.0 ml-collections==0.1.0 numpy==1.21.6 scipy==1.7.0 protobuf==3.20.1 pandas==1.3.4 tensorflow==2.9.0 tensorflow-cpu==2.9.0 matplotlib==3.6.2 python-igraph==0.9.10 pyyaml future csb psutil paramiko scikit-learn pickle5 jinja2
pip install jax==0.3.25 jaxlib==0.3.25+cuda11.cudnn805 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
pip install git+https://github.com/deepmind/alphafold.git@v2.2.4
path=$(python -c 'import site; print(site.getsitepackages()[0])')
cd $path
tmpfile=$(mktemp /tmp/openmm.XXXXXX)
curl https://raw.githubusercontent.com/deepmind/alphafold/v2.2.4/docker/openmm.patch -o $tmpfile
patch -p0 < $tmpfile
rm $tmpfile
wget -q -P ${path}/alphafold/common/ https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt
