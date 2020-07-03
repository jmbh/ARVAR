#!/bin/bash

mkdir "$TMPDIR"/ARVARSim/

cd "$HOME"/ARVARSim

sbatch -a 6-100 Sim_jobs.sh