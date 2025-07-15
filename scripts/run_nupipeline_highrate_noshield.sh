#!/bin/sh
# usage: run_pipe.sh PATH/TO/OBSID

# syntax: run_pipe.sh INDIR
if [ $# != 1 ] ; then
    echo Syntax:  run_nupipeline_highrate.sh PATH_TO_OBSID
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

OUTDIR=$INDIR/event_cl_noshield
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

# From right:
# 0th bit: bad pixels from CALDB bad pixel file. Leave that on.
# 1st bit: Falls in bad pixel from on-board disabled pixel file. Leave this on.
# 2nd bit: Event falls in bad pixel from the usrbadpix file. Leave this on.
# 3rd bit: Neighbor of a bad pixel. Ignore this bit.
# 4th bit: Event is on the detector edge. Ignore this bit.
# 5th bit: Event is in a hot / flickering pixel. Ignore this bit.
# 6th bit: Event is in Neighbor is a hot/flickering pixel. Ignore this bit.
# 7th bit: Event fails depth cut. Leave this on.
# 8th bit: Event fails baseline cut. Leave this on.
# 9th bit: Event fails prior/reset cut (FPMA only). Ignore this.
# 10th bit: Event fails prior cut. Ignore this.
# 11th bit: Event Fails reset cut. Ignore this.
# 12th bit: Event with PI out of range (below 0 or above 4096). Leave this on (for now).
type="(STATUS==b0000xxx00xxxx000)"

echo
echo Running pipeline...

cmd=" nupipeline \
clobber=yes \
indir=$INDIR steminput=$STEMINPUTS \
outdir=$OUTDIR \
entrystage=$ENTRYSTAGE exitstage=$EXITSTAGE \
runsplitsc=yes statusexpr=$type \
pntra=OBJECT cleancols=no pntdec=OBJECT"

echo $cmd > $logfile 2>&1
$cmd >> $logfile 2>&1



