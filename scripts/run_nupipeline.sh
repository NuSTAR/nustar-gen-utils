#!/bin/sh
# usage: run_pipe.sh PATH/TO/OBSID

# syntax: run_pipe.sh INDIR
if [ $# != 1 ] ; then
    echo Syntax:  run_nupipeline_solar.sh PATH_TO_OBSID
    exit 1
fi

# Note, PATH_TO_OBSID is assumed to be relative to the current
# directory:
INDIR=$1


# Set up your local NuSTAR science environment here:
if [ -z "$HEADAS" ]; then
    echo "Need to initialize the HEASOFT first!"
    exit
fi

OUTDIR=$INDIR/event_cl
if [ ! -d $OUTDIR ]; then
#    echo $OUTDIR needs to be produced
    mkdir -m 750 $OUTDIR
fi

# Assume that INDIR will be the complete path, and we only want the last bit
# for the stem inputs:
STEMINPUTS=nu`basename ${1}`

logfile=$OUTDIR/$$_pipe.log

# Set the entry/exit stages here if you want to 
# change it from the default of 1 and 2, respectively.
# Only set EXISTAGE=3 if you actually know what you're doing and have
# added correct keywords for effective area, grprmf, vignetting, etc below.

ENTRYSTAGE=1
EXITSTAGE=2

echo
echo Running pipeline...

cmd=" nupipeline \
clobber=yes \
indir=$INDIR steminput=$STEMINPUTS \
outdir=$OUTDIR \
entrystage=$ENTRYSTAGE exitstage=$EXITSTAGE \
runsplitsc=yes \
pntra=OBJECT cleancols=no pntdec=OBJECT"

echo $cmd > $logfile 2>&1
$cmd >> $logfile 2>&1



