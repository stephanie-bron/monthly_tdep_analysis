#_LAB_MAIN_DIR=$PWD
_LAB_MAIN_DIR=/data/user/sbron/monthly_tdep_analysis_test/
#_LAB_MAIN_DIR='psLab_rev140020_motnly_ana61756a84-4ef3-4841-800a-86a17e866aed'
_LAB_LIB_DIR=$_LAB_MAIN_DIR/psLab/lib64
ICETRAY_ENV=/data/user/achristov/meta-projects/icerec/V04-11-02/build/env-shell.sh
_LAB_COORDINATE_PROJECT=coord_interface

_LAB_DATA_LOCATION=Geneve

_LAB_DATA_DIR=/data/user/achristov/data
_LAB_RESULTS_DIR=/data/user/achristov/results

_LAB_CORE_DIR=$_LAB_MAIN_DIR/psLab/labcore
_LAB_PRELOAD=$_LAB_CORE_DIR/preload.C

if [[ -z "$LAB_MAIN_DIR" ]]
    then
    
    printf "\nWelcome to the lab.\n"
    printf "   LAB_MAIN_DIR = %s\n" $_LAB_MAIN_DIR
    printf "   LAB_LIB_DIR  = %s\n" $_LAB_LIB_DIR
    printf "   LAB_CORE_DIR = %s\n" $_LAB_CORE_DIR
    printf "   LAB_PRELOAD  = %s\n" $_LAB_PRELOAD
    printf "   LAB_DATA_LOCATION = %s\n" $_LAB_DATA_LOCATION
    printf "   LAB_DATA_DIR = %s\n" $_LAB_DATA_DIR
    printf "   ICETRAY_ENV  = %s\n" $ICETRAY_ENV
    printf "   LAB_COORDINATE_PROJECT = %s\n" $_LAB_COORDINATE_PROJECT
    printf "\n"

    LAB_MAIN_DIR=$_LAB_MAIN_DIR \
	LAB_LIB_DIR=$_LAB_LIB_DIR \
	LAB_CORE_DIR=$_LAB_CORE_DIR \
	LAB_PRELOAD=$_LAB_PRELOAD \
	LAB_DATA_LOCATION=$_LAB_DATA_LOCATION \
        LAB_DATA_DIR=$_LAB_DATA_DIR \
	LAB_COORDINATE_PROJECT=$_LAB_COORDINATE_PROJECT \
	$ICETRAY_ENV $@ 
          # $@ allows passing arguments to the new shell

    printf "You are exiting the lab.\n"
    exit
else
    printf "You are already in the lab.  You must exit before reloading.\n"
fi
