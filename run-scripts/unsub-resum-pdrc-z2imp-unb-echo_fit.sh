#!/bin/bash -l

SBATCH_FILE="./fits/echo_sbatch_scripts/uns-res-pdrc-z2imp-unb-FIT_Qs0_C2_-cueps-$1-cumeval-$2-minuprec-$3-fitpar-$4.sh"
echo $SBATCH_FILE

INPUT_FILE="./fits/echo_inputs/unsub-resum-pdrc-z2imp-unb-fit.sh"

echo "#!/bin/bash -l" > $SBATCH_FILE
echo "#SBATCH -J f2-URPI # unsub resum parent (z2)improved" >> $SBATCH_FILE
echo "#SBATCH -o ./fits/out/unsub-resumbk-pdrc-z2imp-unboundloop/fit_NUM_TESTING_FIT_Qs0_C2_cueps-$1-cumeval-$2-minuprec-$3-fitpar-$4_%j.txt" >> $SBATCH_FILE
cat $INPUT_FILE >> $SBATCH_FILE

echo "# fit NUM settings" >> $SBATCH_FILE
CUBAMTHD="suave"
echo "CUBAMTHD="'"'"$CUBAMTHD"'"' >> $SBATCH_FILE
EPS=$1
echo "EPS=$EPS" >> $SBATCH_FILE
MEVAL=$2
echo "MEVAL=$MEVAL" >> $SBATCH_FILE
MINUITPREC=$3
echo "MINUITPREC=$MINUITPREC" >> $SBATCH_FILE
FITPAR=$4
echo "FITPAR=$FITPAR" >> $SBATCH_FILE

echo 'echo FIT HERA DATA $SCHEME $BK $RC $Z2BOUND $BOUNDLOOP $QSSQR $CSQR $GAMMA $X0IF $X0BK $EC $TYPVIRTQ02 $Y0 $ETA0 $CUBAMTHD $EPS $MEVAL $MINUITPREC $FITPAR' >> $SBATCH_FILE
echo 'CUBACORES=0 ./fitex $SCHEME $BK $RC $Z2BOUND $BOUNDLOOP $QSSQR $CSQR $GAMMA $X0IF $X0BK $EC $TYPVIRTQ02 $Y0 $ETA0 $CUBAMTHD $EPS $MEVAL $MINUITPREC $FITPAR' >> $SBATCH_FILE
#echo 'seff $SLURM_JOBID' >> $SBATCH_FILE

#cat $SBATCH_FILE
#actual final command
sbatch $SBATCH_FILE


