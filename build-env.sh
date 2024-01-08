#!/bin/bash

# build-env.sh copied from:
# https://repo.maap-project.org/sshah/gdal-hello-world/-/blob/main/build-env.sh

# source activate s1rtc
basedir=$( cd "$(dirname "$0")" ; pwd -P )
# conda install ${basedir}/environment.yaml
conda config --set solver libmamba

# conda env update --file ${basedir}/environment.yaml
conda env create -f environment.yaml
source activate s1rtc
conda install -y -c conda-forge numpy=1.25
