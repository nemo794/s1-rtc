#!/bin/bash
# This is intended for running DPS jobs.
# For each input where download=TRUE in the algorithm_config.yaml file, 
# DPS places it into the "input" directory for the algorithm to use.

# Example run.sh scripts:
# https://repo.maap-project.org/sshah/gdal-hello-world/-/blob/main/run.sh
# https://github.com/lauraduncanson/icesat2_boreal/blob/08e36c172a604013c6f69f6b917379dbf4c893c3/dps/alg_3-1-5/run.sh

# For custom environments, make sure to activate the conda environment.
# For DPS, use `source activate <custom>`, not `conda activate <custom>`
source activate s1rtc

# TODO - move this section to a build-env.sh script
# basedir=$( cd "$(dirname "$0")" ; pwd -P )
# # conda install ${basedir}/environment.yaml
# mamba env update --file ${basedir}/environment.yaml

# INPUT_FILENAME=$(ls -d input/*)
INPUT_FILENAME=$1

# Get path to this run.sh script
basedir=$( cd "$(dirname "$0")" ; pwd -P )

# Per DPS convention, the directory to place outputs into MUST be called "output".
# Only items in a directory with that name will persist in my-public-bucket after DPS finishes.
mkdir -p output

export HOME=/home/ops

python ${basedir}/s1_rtc.py --in_file ${INPUT_FILENAME} --output_dir output
# python ${basedir}/s1_rtc.py --in_file input/${INPUT_FILENAME} --output_dir output
