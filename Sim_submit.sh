#!/bin/bash

mkdir "$TMPDIR"/ARVARSim/

cd "$HOME"/ARVARSim

for i in `seq 1 100`;
do
	sed s/iter/$i/g Sim_jobs.sh > CUR_submit.sh
    sbatch CUR_submit.sh
done

rm -r "$TMPDIR"/ARVARSim/
rm -r CUR_submit.sh



