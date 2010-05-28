#!/bin/bash
LEVEL='DEBUG'
#LEVEL='INFO'

nevt=$1
if [ -z $nevt ]
then
  nevt=1
fi

time rivet -a FASTSIM -H fastsim -x 1 -n ${nevt}  -l Rivet.Analysis=${LEVEL} events.hepmc


