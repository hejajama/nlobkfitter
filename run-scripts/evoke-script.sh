#!/bin/bash -l

#SettingList="unsub-resum-pdrc-z2imp-unb unsub-kcbk-pdrc-z2imp-unb unsub-resum-sdrc-z2imp-unb unsub-kcbk-sdrc-z2imp-unb unsub-trbk-pdrc-z2sim-unb unsub-trbk-sdrc-z2sim-unb"
SettingList="unsub-resum-pdrc-z2imp-unb unsub-resum-sdrc-z2imp-unb"
#SettingList="unsub-kcbk-pdrc-z2imp-unb unsub-kcbk-sdrc-z2imp-unb"
#SettingList="unsub-trbk-pdrc-z2sim-unb unsub-trbk-sdrc-z2sim-unb"

# mkdir loop
#for name in $SettingList; do
#  mkdir ./out/$name
#done
#exit 1

for SETTING in $SettingList; do
  for cueps in 1e-3 1e-4; do # 
    for cumeval in 2e7 4e7; do # 
      for minuprec in 1e-4; do # 
        for fitpar in "01" "2"; do
          echo model: $SETTING eps: $cueps , cumeval: $cumeval, minuprec: $minuprec, fitpar: $fitpar
          sh ./fits-lightq/echo-fit-sbatch.sh $SETTING $cueps $cumeval $minuprec $fitpar
#          echo ./fits-lightq/echo-fit-sbatch.sh $SETTING $cueps $cumeval $minuprec $fitpar
        done
      done
    done
  done
done
