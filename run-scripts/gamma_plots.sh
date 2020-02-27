#!/bin/bash

# COMPUTE DIPOLES AND STRUCTURE FUNCTIONS IN LHeC KINEMATICS
# Compare pdrc and sdrc for KCBK

# HERA FITS

# x0bk = 0.01
    # pdrc - final fit
./plottool unsub kcbk pdrc z2imp unb 0.0833 3.49326 0.98 0.01 25.0109 | tee ./out/lhec-fin-hera-ukpi-x001.dat
    # sdrc - final fit
./plottool unsub kcbk sdrc z2imp unb 0.0905 0.845969 1.21 0.01 22.2806 | tee ./out/lhec-fin-hera-uksi-x001.dat

# x0bk = 1.0
    # pdrc - temp
#./plottool unsub kcbk pdrc z2imp unb 0.085 4.24795 1.0 0.01 26.5772 | tee ./out/lhec-hera-ukpi-x1.dat
    # sdrc - temp
#./plottool unsub kcbk sdrc z2imp unb 0.085 4.24795 1.0 0.01 26.5772 | tee ./out/lhec-hera-uksi-x1.dat

