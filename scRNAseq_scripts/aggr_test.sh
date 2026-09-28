#!/bin/bash
#SBATCH -J aggr
#SBATCH -o aggr.out
#SBATCH -N 1 -n 32

module load cellranger

echo "sample_id,molecule_h5" > control_aggr.csv
echo "control_1,./control_1/outs/molecule_info.h5" >> control_aggr.csv
echo "control_2,./control_2/outs/molecule_info.h5" >> control_aggr.csv

cellranger aggr --id=control_aggr --csv=control_aggr.csv

echo "sample_id,molecule_h5" > treatment_aggr.csv
echo "treatment_1,./treatment_1/outs/molecule_info.h5" >> treatment_aggr.csv
echo "treatment_2,./treatment_2/outs/molecule_info.h5" >> treatment_aggr.csv

cellranger aggr --id=treatment_aggr --csv=treatment_aggr.csv

