#!/bin/bash -i

echo 'ARCIMBOLDO_AIR INSTALLATION'
echo 'The script will proceed to verify that all the dependencies have been installed and install the remaining ones'
echo 'Checking conda...'
if command -v conda &> /dev/null; then
    echo "Conda is installed."
    if conda info &> /dev/null; then
        echo "Conda is executable."
    else
        echo "Conda is installed but cannot be executed. Verify Conda can be executed and relaunch the installation script."
        exit 1
    fi
else
    echo "Conda is not installed (https://conda.io/projects/conda/en/stable/user-guide/install/download.html)."
    # Ask user if they want to install Conda
    read -p "Do you want to install Conda? (y/n): " choice
    if [ "$choice" == "y" ] || [ "$choice" == "Y" ]; then
        curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
        bash Miniconda3-latest-Linux-x86_64.sh
        rm Miniconda3-latest-Linux-x86_64.sh
        echo "Conda has been installed. Please open a new terminal session to use Conda."
    else
        echo "Conda installation skipped. Exiting the installation script."
    fi
    exit 1
fi

# Check if Nvidia driver module is loaded
echo 'Checking Nvidia driver...'
if lsmod | grep -i nvidia &> /dev/null; then
    echo "Nvidia driver is installed and loaded."
else
    echo "Nvidia driver is not installed or loaded. Please, install Nvidia driver before continuing (https://www.nvidia.es/Download/index.aspx?lang=eng)."
    #exit 1
fi

# Check if maxit command exists
if command -v maxit &> /dev/null; then
    echo "Maxit is installed."
else
    echo "Maxit is not installed. In order to continue, download and install maxit (https://sw-tools.rcsb.org/apps/MAXIT/index.html)."
    exit 1
fi

# Check CCP4
if command -v pisa &> /dev/null && command -v pdbset &> /dev/null && command -v lsqkab &> /dev/null; then
    echo "CCP4 suite is installed."
else
    echo "CCP4 programs could not be found in the PATH. Please, check if lsqkab, pisa and pdbset can be found in the PATH."
fi

echo 'Creating ARCIMBOLDO_AIR Conda environment.'
read -p "Enter the Conda environment name: " env_name
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda create -y -n "$env_name" python=3.8
conda activate "$env_name"
conda install -y -c conda-forge openmm==7.5.1 cudatoolkit==11.2 cudnn==8.2.1.32 pdbfixer==1.7
conda install -y -c bioconda hmmer==3.3.2 hhsuite==3.3.0 kalign2==2.04

# Check if ptxas command exists
echo 'Checking ptxas...'
if command -v ptxas &> /dev/null; then
    echo "ptxas command is installed. Skipping cudatoolkit-dev installation"
else
  if conda install -y -c conda-forge "cudatoolkit-dev=11.2"; then
      echo 'Cudatoolkit-dev installed.'
  else
      echo 'Conda installation of cudatoolkit failed. Please, install cudatoolkit-dev afterwards.'
  fi
fi

pip install absl-py==1.0.0 biopython==1.79 chex==0.0.7 dm-haiku==0.0.9 dm-tree==0.1.6 immutabledict==2.0.0 ml-collections==0.1.0 numpy==1.21.6 scipy==1.7.0 protobuf==3.20.1 pandas==1.3.4 tensorflow==2.9.0 tensorflow-cpu==2.9.0 matplotlib==3.6.2 python-igraph==0.9.10 pyyaml future csb psutil paramiko scikit-learn pickle5 jinja2
pip install jax==0.3.25 jaxlib==0.3.25+cuda11.cudnn805 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
pip install git+https://github.com/deepmind/alphafold.git@v2.2.4
path=$(python -c 'import site; print(site.getsitepackages()[0])')
cd "$path" || exit
tmpfile=$(mktemp /tmp/openmm.XXXXXX)
curl https://raw.githubusercontent.com/deepmind/alphafold/v2.2.4/docker/openmm.patch -o "$tmpfile"
patch -p0 < "$tmpfile"
rm "$tmpfile"
wget -q -P "${path}"/alphafold/common/ https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt
