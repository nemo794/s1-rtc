#!/bin/bash

# build-env.sh copied from:
# https://repo.maap-project.org/sshah/gdal-hello-world/-/blob/main/build-env.sh

source activate base
basedir=$( cd "$(dirname "$0")" ; pwd -P )
conda install ${basedir}/environment.yaml
