#SBATCH -J f2-URSI # unsub resum smallest (z2)improved
#SBATCH -e /dev/null
##SBATCH -n 8
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=90000
##SBATCH -t 12:00:00
#SBATCH -t 72:00:00
##SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=hejohann@student.jyu.fi
#

SCHEME=unsub            # sub / unsub / unsub+
BK=resumbk              # resumbk / kcbk / lobk
RC=smallestrc             # parentrc / guillaumerc / fixedrc
Z2BOUND=z2improved        # z2improved / z2simple
BOUNDLOOP=unboundloop   # z2boundloop / unboundloop

# fit IC
QSSQR=0.05
CSQR=1.36
GAMMA=2.8
X0IF=1.0
X0BK=1.0
EC=1.0
TYPVIRTQ02=1.0
Y0=0
ETA0=0

