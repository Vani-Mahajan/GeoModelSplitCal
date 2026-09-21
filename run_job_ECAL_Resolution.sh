#!/bin/bash

PROCESS=$1

# Gap array




ENERGIES=(
5000 10000 20000 50000 100000
)



ENERGY=${ENERGIES[$PROCESS]}


source /cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc13-opt/setup.sh

SCRATCH=$_CONDOR_SCRATCH_DIR

mkdir -p $SCRATCH

cp -r /afs/cern.ch/user/m/mahajanv/GeoModelSplitCal $SCRATCH/

cd $SCRATCH/GeoModelSplitCal


RUN_FILE="run.cfg"
sed -i "s/^.*energy_MeV .*$/energy_MeV = ${ENERGY}/" ${RUN_FILE}

# Copy macro and config files into build directory so Geant4 can access them
cp run.cfg build/
cp run.mac build/

cd $SCRATCH/GeoModelSplitCal/build

chmod +x run_g4

./run_g4

OUT_DIR="/afs/user/m/mahajanv/GeoModelSplitCal/sim_results"
mkdir -p $OUT_DIR

# Safely find and move the output ROOT file
ROOT_FILE=$(find . -maxdepth 1 -type f -name "*.root" | head -n 1)

if [ -n "$ROOT_FILE" ]; then
    echo "Found ROOT file: $ROOT_FILE"
    mv "$ROOT_FILE" "${OUT_DIR}/sim_e-_10k_${ENERGY}MeV.root"
    echo "Saved:"
    echo "${OUT_DIR}/sim_e-_10k_${ENERGY}MeV.root"
else
    echo "ERROR: No ROOT file produced!"
    exit 1
fi


