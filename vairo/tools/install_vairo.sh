#!/bin/bash -i

echo "VAIRO INSTALLATION"
echo "The script will verify the installation of all dependencies and proceed to install any remaining ones."

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
        echo "Conda has been installed. Please open a new terminal session and rerun the script."
    else
        echo "Conda installation skipped. Exiting the installation script."
    fi
    exit 1
fi

# Check if Nvidia driver module is loaded
echo "Checking Nvidia driver..."
if lsmod | grep -i nvidia &> /dev/null; then
    echo "Nvidia driver is installed and loaded."
else
    echo "Nvidia driver is not installed or loaded. Please, install Nvidia driver before continuing (https://www.nvidia.es/Download/index.aspx?lang=eng)."
    exit 1
fi

# Check if maxit command exists
echo 'Checking MAXIT...'
if command -v maxit &> /dev/null; then
    echo "MAXIT is installed."
else
    echo "MAXIT is not installed. In order to continue, download and install maxit (https://sw-tools.rcsb.org/apps/MAXIT/index.html)."
    exit 1
fi

# Check CCP4
echo 'Checking CCP4...'
if command -v pisa &> /dev/null && command -v pdbset &> /dev/null && command -v lsqkab &> /dev/null; then
    echo "CCP4 suite is installed."
else
    echo "The CCP4 programs cannot be found in the PATH. Please check if LSQKAB, PISA and PDBSET are present in the PATH. If the CCP4 suite is not installed, please install it from the following link: https://www.ccp4.ac.uk/"
fi

echo 'Creating VAIRO Conda environment.'
read -p "Enter the Conda environment name: " env_name
conda create -y -n "$env_name" python=3.8
conda activate "$env_name"
conda install -y -c conda-forge openmm==7.5.1 cudatoolkit==11.2 cudnn==8.2.1.32 pdbfixer==1.7
conda install -y -c bioconda hmmer==3.3.2 hhsuite==3.3.0 kalign2==2.04

# Check if ptxas command exists
echo "Checking ptxas..."
if command -v ptxas &> /dev/null; then
    echo "ptxas command is installed. Skipping cudatoolkit-dev installation."
else
  if conda install -y --quiet -c conda-forge "cudatoolkit-dev=11.2"; then
      echo "Cudatoolkit-dev installed."
  else
      echo "Conda installation of cudatoolkit-dev failed. In order to continue, download and install cudatoolkit: https://developer.nvidia.com/cuda-toolkit"
      exit 1
  fi
fi

pip install absl-py==1.0.0 biopython==1.79 chex==0.0.7 dm-haiku==0.0.9 dm-tree==0.1.6 immutabledict==2.0.0 ml-collections==0.1.0 numpy==1.21.6 scipy==1.7.0 protobuf==3.20.1 pandas==1.3.4 tensorflow==2.9.0 tensorflow-cpu==2.9.0 matplotlib==3.6.2 python-igraph==0.9.10 pyyaml future csb psutil paramiko scikit-learn pickle5 jinja2 flask
pip install jax==0.3.25 jaxlib==0.3.25+cuda11.cudnn805 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
pip install git+https://github.com/deepmind/alphafold.git@v2.2.4
path=$(python -c "import site; print(site.getsitepackages()[0])")
cd "$path" || exit
tmpfile=$(mktemp /tmp/openmm.XXXXXX)
curl https://raw.githubusercontent.com/deepmind/alphafold/v2.2.4/docker/openmm.patch -o "$tmpfile"
patch -p0 < "$tmpfile"
rm "$tmpfile"
wget -q -P "${path}"/alphafold/common/ https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

echo "******"
echo "VAIRO HAS BEEN SUCCESSFULLY INSTALLED"
echo "A conda environment with the name ""$env_name"" has been created. To run VAIRO, you must first activate the new environment (conda activate ""$env_name"") and type VAIROGUI. This will run the program and show the different options for creating a job."
echo "Before running VAIRO, please ensure that the libraries of AlphaFold2 are installed. If they have not been installed, please download them using the following script.: https://github.com/deepmind/alphafold/blob/v2.2.4/scripts/download_all_data.sh"
echo "Run the following commands:"
echo "conda activate ""$env_name"""
echo "VAIROGUI"