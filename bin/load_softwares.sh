#!/bin/bash

# Load environment variables
source project_env.sh

# Define and create the software directory
cd "${my_softwares}"

# Clone SEACR if it doesn't already exist
if [ ! -d SEACR ]; then
  echo "📥 Cloning SEACR..."
  git clone https://github.com/FredHutch/SEACR.git
  chmod +x SEACR/SEACR_1.3.sh
else
  echo "✅ SEACR already installed."
fi

# Ensure SEACR script is executable
chmod +x SEACR/SEACR_1.3.sh

module load conda
conda_env_path="${my_softwares}/conda_envs/bigwig"
mkdir -p "$(dirname "$conda_env_path")"

if [ ! -d "$conda_env_path" ]; then
  echo "📦 Creating Conda env 'bigwig' under ${conda_env_path}..."
  conda create --prefix "$conda_env_path" ucsc-bigwigtobedgraph -c bioconda -y
else
  echo "✅ Conda environment already exists at ${conda_env_path}"
fi

# Step 5: Activate and test
echo "To activate:"
echo "  conda activate $conda_env_path"