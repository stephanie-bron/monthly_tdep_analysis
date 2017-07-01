#!/bin/zsh
export HEADAS=/data/user/achristov/FERMI/tools/heasoft/x86_64-unknown-linux-gnu-libc2.12/
#export HEADAS=/data/user/sbron/FERMI/tools/heasoft/x86_64-unknown-linux-gnu-libc2.12/
. $HEADAS/headas-init.sh

export FERMI_DIR=/data/user/achristov/FERMI/tools/ScienceTools-v10r0p5-fssc-20150518-i686-pc-linux-gnu-libc2.5-without-root/i686-pc-linux-gnu-libc2.5/
#export FERMI_DIR=/data/user/sbron/FERMI/tools/ScienceTools-v10r0p5-fssc-20150518-i686-pc-linux-gnu-libc2.5-without-root/i686-pc-linux-gnu-libc2.5/
source $FERMI_DIR/fermi-init.sh
eval `/cvmfs/icecube.opensciencegrid.org/py2-v1/setup.sh`
export PYTHONPATH=$PYTHONPATH:/data/user/achristov/monthly_tdep/astropy

#astropy installed using
#pip install --target=/data/user/achristov/monthly_tdep/astropy astropy

#alabaster installed using
#pip install --target=/data/user/achristov/monthly_tdep/alabaster  alabaster