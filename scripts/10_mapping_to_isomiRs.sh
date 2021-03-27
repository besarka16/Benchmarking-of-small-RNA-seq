#!/usr/bin/env bash

if [ $# -eq 0 ]
  then
    echo "Supply number of cores!"
fi

CORES=$1


working_dir=$(pwd)

cd isomiRROR-miRXplore/
sudo -E env "PATH=$PATH" snakemake --use-conda --cores $CORES

cd $working_dir

cd isomiRROR-plasma
sudo -E env "PATH=$PATH" snakemake --use-conda --cores $CORES

cd $working_dir

cp isomiRROR-plasma/isomir_readcount.txt results/plasma_isomir_readcount.txt
cp isomiRROR-miRXplore/isomir_readcount.txt results/miRXplore_isomir_readcount.txt