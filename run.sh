#!/bin/bash
# This is intended for running DPS jobs.
# For each input where download=TRUE in the algorithm_config.yaml file, 
# DPS places it into the "input" directory for the algorithm to use.

# Example run.sh scripts:
# https://repo.maap-project.org/sshah/gdal-hello-world/-/blob/main/run.sh
# https://github.com/lauraduncanson/icesat2_boreal/blob/08e36c172a604013c6f69f6b917379dbf4c893c3/dps/alg_3-1-5/run.sh

# For custom environments, make sure to activate the conda environment.
# For DPS, use `source activate <custom>`, not `conda activate <custom>`
source activate base

INPUT_FILENAME=$(ls -d input/*)
# INPUT_FILENAME=$1

# Get path to this run.sh script
basedir=$( cd "$(dirname "$0")" ; pwd -P )

# Per DPS convention, the directory to place outputs into MUST be called "output".
# Only items in a directory with that name will persist in my-public-bucket after DPS finishes.
mkdir -p output

python ${basedir}/s1_rtc.py --in_file ${INPUT_FILENAME} --output_dir output
# python ${basedir}/s1_rtc.py --in_file input/${INPUT_FILENAME} --output_dir output

# set -x
# unset PROJ_LIB

# # https://stackoverflow.com/questions/42352841/how-to-update-an-existing-conda-environment-with-a-yml-file
# # https://stackoverflow.com/questions/36539623/how-do-i-find-the-name-of-the-conda-environment-in-which-my-code-is-running
# conda env update --name $CONDA_DEFAULT_ENV --file environment.yaml

# # Absolute path here
# # This PWD is wherever the job is run (where the .sh is called from) 
# OUTPUTDIR="${PWD}/output"

# # Cmd line call that worked
# #python 3.1.5_dps.py --in_tile_fn '/projects/maap-users/alexdevseed/boreal_tiles.gpkg' --in_tile_num 30550 --tile_buffer_m 120 --in_tile_layer "boreal_tiles_albers" -o '/projects/tmp/Topo/'

# #Print to stdout for debugging
# python s1_rtc.py \
# --s1_zip $INPUT_FILENAME \
# --output_dir $OUTPUTDIR \
