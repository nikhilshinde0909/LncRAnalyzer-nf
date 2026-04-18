#!/bin/bash

# Check for existing Mambaforge, Miniforge, or Anaconda installation
if [[ -d "$HOME/mambaforge" ]]; then
    echo "Existing Mambaforge installation detected."
    BIN=$HOME/mambaforge/bin/
elif [[ -d "$HOME/miniforge" ]]; then
    echo "Existing Miniforge installation detected."
    BIN=$HOME/miniforge/bin/
elif [[ -d "$HOME/anaconda3" ]]; then
    echo "Existing Anaconda installation detected."
    BIN=$HOME/anaconda3/bin/
else
    echo "No recognized environment (Mambaforge, Miniforge, Anaconda) found in $HOME."
    echo "Installing Miniforge..."

    curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -o miniforge.sh \
    && chmod +x miniforge.sh \
    && bash miniforge.sh -b -p $HOME/miniforge \
    && rm miniforge.sh

    BIN=$HOME/miniforge/bin/
fi

# Export paths
export PATH=$PATH:$BIN

# Verify that mamba is installed, if not use conda
if ! command -v mamba &> /dev/null; then
    echo "Mamba not found, trying conda..."
    mamba="conda"
else
    mamba="mamba"
fi

# Update the mamba environment using the YAML file
echo "Updating existing environment..."
$mamba env update --file LncRAnalyzer.yml && conda clean -y -a

# Installing LncFinder
if [[ -f "$BIN/R" ]]; then
    echo "R installation found.."
    R="$BIN/R"
    $R -e 'options(timeout=600); install.packages(c("LncFinder", "seqinr"), repos="https://cran.r-project.org")'
else
    echo "R installation not found..."
    exit 1
fi

# Creating FEELnc environment
echo "Creating FEELnc environment..."
$mamba create -y -n FEELnc -c bioconda feelnc && conda clean -y -a

# Creating CPC-CPAT-Slncky environment
echo "Creating CPC-CPAT-Slncky environment..."
$mamba env create -f cpc2-cpat-slncky.yml && conda clean -y -a

# Creating RNAsamba environment
echo "Creating RNAsamba environment..."
$mamba env create -f rnasamba.yml && conda clean -y -a

# install hmmer=3.1b1 from source
echo "Installing hmmer=3.1b1 from source"

function hmmer_install {
    curl -L -o hmmer-3.1b1.tar.gz http://eddylab.org/software/hmmer/hmmer-3.1b1.tar.gz
    tar zxvf hmmer-3.1b1.tar.gz ; rm -rf hmmer-3.1b1.tar.gz
    mv hmmer-3.1b1 utils/
    cd utils/hmmer-3.1b1/
    ./configure
    make
    cd ..
    cd ..
}

hmmer_install

# Create the symbolic link for hmmmer-3.1b
ln -sf $PWD/utils/hmmer-3.1b1/src/* $BIN
hmmer_path=`which hmmscan 2>/dev/null`
echo "HMMER 3.1b1 has been installed to: ${hmmer_path}"

# Detect existing Python2.7 installation and change permissions
if [[ -d "$HOME/mambaforge/envs/cpc2-cpat-slncky/" ]]; then
    echo "Python2.7 installation detected in Mambaforge."
    PY_DIR="$HOME/mambaforge/envs/cpc2-cpat-slncky/lib/python2.7"
    LINK_PATH="$HOME/mambaforge/envs/cpc2-cpat-slncky/bin/"
elif [[ -d "$HOME/miniforge/envs/cpc2-cpat-slncky/" ]]; then
    echo "Python2.7 installation detected in Miniforge."
    PY_DIR="$HOME/miniforge/envs/cpc2-cpat-slncky/lib/python2.7"
    LINK_PATH="$HOME/miniforge/envs/cpc2-cpat-slncky/bin/"
elif [[ -d "$HOME/anaconda3/envs/cpc2-cpat-slncky/" ]]; then
    echo "Python2.7 installation detected in Anaconda."
    PY_DIR="$HOME/anaconda3/envs/cpc2-cpat-slncky/lib/python2.7"
    LINK_PATH="$HOME/anaconda3/envs/cpc2-cpat-slncky/bin/"
else
    echo "No Python2.7 found in $HOME."
    exit 1
fi
echo "Changing permissions for $PY_DIR."
chmod -R 777 $PY_DIR
echo "Creating symbolic link for Slncky"
chmod +x $PWD/utils/slncky/slncky.v1.0 $PWD/utils/slncky/alignTranscripts1.0
ln -sf $PWD/utils/slncky/alignTranscripts1.0 $PWD/utils/slncky/slncky.v1.0 $LINK_PATH

echo "Installation complete !!"
