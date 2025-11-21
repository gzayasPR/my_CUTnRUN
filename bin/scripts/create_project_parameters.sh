#!/bin/bash

# Script: create_project_parameters.sh
# Purpose: Generate a standardized project-wide parameter file from user input
# Usage:   ./create_project_parameters.sh Thermo_Seminole.parameter.sh

set -euo pipefail

INPUT_PARAM_FILE="$1"

# --- 1. Load input parameter file ---
if [[ ! -f "$INPUT_PARAM_FILE" ]]; then
    echo "ERROR: Input parameter file '$INPUT_PARAM_FILE' not found."
    exit 1
fi

# Load user-defined variables
source "$INPUT_PARAM_FILE"

# --- 2. Check required variables (edit as needed) ---
required_vars=(
    NAME
    CSV
    account
    ECOLI_INDEX
    MOUSE_INDEX
    MOUSE_GTF
    MOUSE_BED
    BLACKLIST
    my_results
    my_bin
    my_data
    my_bash
)

for var in "${required_vars[@]}"; do
    if [[ -z "${!var:-}" ]]; then
        echo "ERROR: Required variable '$var' is not set in $INPUT_PARAM_FILE"
        exit 2
    fi
done

# --- 2b. Set default if missing ---
: "${ONE_LV_PER_GROUP:=true}"

# --- 3. Set up output directory for parameters file ---
PARAM_DIR="${my_results}/${NAME}/parameters"
mkdir -p "$PARAM_DIR"
PARAM_FILE="${PARAM_DIR}/project_env.parameters.sh"

# --- 4. Write standardized parameter file ---
cat << EOF > "$PARAM_FILE"
#############################################################
#  Auto-generated project parameter file for: $NAME
#  Generated from $INPUT_PARAM_FILE
#  Created: $(date)
#############################################################

# -- Project & File Paths --
export NAME="${NAME}"
export CSV="${CSV}"
export account="${account}"
export BLACKLIST="${BLACKLIST}"
export ECOLI_INDEX="${ECOLI_INDEX}"
export MOUSE_INDEX="${MOUSE_INDEX}"
export MOUSE_GTF="${MOUSE_GTF}"
export MOUSE_BED="${MOUSE_BED}"
export OUT_DIR="${my_results}/${NAME}"
export bash_out="${my_bash}/${NAME}"

export my_bin="${my_bin}"
export my_data="${my_data}"
export my_results="${my_results}"
export my_bash="${my_bash}"

# -- SLURM Resource Parameters --
export STEP1_CPUS=\${STEP1_CPUS:-1}
export STEP2_CPUS=\${STEP2_CPUS:-2}
export STEP3_CPUS=\${STEP3_CPUS:-4}
export STEP4_CPUS=\${STEP4_CPUS:-4}
export STEP5_CPUS=\${STEP5_CPUS:-4}


#############################################################
EOF

echo "✅ Created parameter file: $PARAM_FILE"
