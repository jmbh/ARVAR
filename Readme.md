# Reproducibility Archive
This Readme describes the files in the reproducibility archive for Dablander, F.<sup><small>&#11089;</small></sup>, Ryan, O.<sup><small>&#11089;</small></sup>, & Haslbeck, J.M.B.<sup><small>&#11089;</small></sup> (under review). Choosing between AR(1) and VAR(1) Models in Typical Psychological Applications. [[Link](https://psyarxiv.com/qgewy/)]. It allows one to reproduce all analyses and figures in the paper.

- SampleGrid.R
    - Input: Mixed Model estimated on "MindMaastricht" data (in /files)
    - Output: 74 x 100 models, 100 in each "cell"
- aux_funcions.R
    - Contains auxilliary functions used for the simulation script and to process the simulated data
- aux_funcions.R
    - Contains auxilliary functions used to visualize the simulated data
- Sim_LISA.R
    - This file runs a single iteration (of 100) of the simulation reported in the paper; it imports functions from aux_functions.R
    - Input: 7400 VAR models in file Models_60cells.RDS
    - Output: 1 RDS output file
- Sim_jobs.sh and Sim_submit.sh
    - These two bash files were used to run Sim_LISA.R 100 times in parallel on the LISA cluster computer at the University of Amsterdam
- Analysis.Rmd
    - Analysis script walking through the creation of all Figures in the paper
- Analysis.html
    - Compiled from Analysis.Rmd

The scripts above make use of the following files:

- MM_Bringmann.RDS: contains the estimated mixed model from the "MindMaastricht" data
- initial10000.RDS: contains R/D values of the initial 10000 VAR models sampled in SampleGrid.R
- Models_74cells.RDS: A list with 74 entries, each of which contains 100 VAR matrices; these are the VAR matrices we use in the simulation file Sim_LISA.R
- cell_positions.RDS: A matrix that contains the cell positions (x/y positions of borders)
- EE_comp.RDS: Finished computation of EE_comp values; included since it takes a while to compute it; this file is created in PlottingFigures.R
- out_var.RDS: List holding all simulation results
- diffs_mat.RDS: A matrix storing the differences in estimation error between AR and VAR

We did not include the 100 simulation output files in this archive due to their size (~15GB). However, we are happy to provide those files on request.


