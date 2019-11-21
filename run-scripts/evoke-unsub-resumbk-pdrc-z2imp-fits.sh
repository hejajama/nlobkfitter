#!/bin/bash -l
  
# run1
#QDN=50
#QUP=300
#QSTEP=50
#CDN=100
#CUP=400
#CSTEP=50



#for cueps in 1e-3 1e-4 1e-5 1e-6; do # run1 Qs0
for cueps in 1e-4 1e-6; do # run2 Qs0 C^2
  #for cumeval in 2e7 4e7 6e7; do # run1
  for cumeval in 2e7 4e7; do # run2
    #for minuprec in 1e-4 1e-5 1e-6 1e-7; do # run1
    for minuprec in 1e-4 1e-5; do # run2
      for fitpar in "01" "1"; do
        echo eps: $cueps , cumeval: $cumeval, minuprec: $minuprec, fitpar: $fitpar
        sh ./fits/unsub-resum-pdrc-z2imp-unb-echo_fit.sh $cueps $cumeval $minuprec $fitpar
      done
    done
  done
done

