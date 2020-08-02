#!/bin/bash
# 
# script to generate and submit independent jobs to HPC cluster.
#
# Author: David Warne (david.warne@qut.edu.au)
#         School of Mathematical Sciences
#         Science and Engineering Faculty
#         Queensland University of Technology
#
# NOTE: must be modified for site specific HPC configurations
#
# NOTE: submits N independent single code PBS jobs, if smc_abc_rw_par is used
#       then set ncpus to a larger number.
#
# NOTE: Due to significant variability in runtimes between countries, the 
#       PBS array job option is not used.

for data in '30mar' '13apr'
do
    outdir="results_$data"
# for every country (with weekly lags for back-testing)
for t in 0
do
    for model in 0
    do
        for (( j=1;j<=250;j++ ))
        do
            # job name
            NAME="T$t""M$model""R$j"
            # checkpoint filename for job
            CHCKFILE="T$t-M$model-R$j-$data-checkpoint.mat"
            
            # write PBS submit script
            cat > sub << EOF
#!/bin/bash -l
#PBS -N $NAME
#PBS -l ncpus=1
#PBS -l walltime=8:00:00
#PBS -l mem=4gb

module load matlab
cd \$PBS_O_WORKDIR
mkdir -p $outdir/T$t/M$model
matlab -r "country_id=$j;data_snap='$data'; pred_horizon=$t;model=$model; OUT_DIR='$outdir/T$t/M$model';checkpointfile='$CHCKFILE'; run_smc_intensity_reg_cluster; exit();"
EOF
            # submit
            qsub sub
            done
        done
    done
done
