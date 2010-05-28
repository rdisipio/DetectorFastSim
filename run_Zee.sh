#!/bin/bash
LEVEL='DEBUG'
#LEVEL='INFO'

nevt=$1
if [ -z $nevt ]
then
  nevt=1
fi

mkfifo ./events.fifo

agile-runmc Pythia6:421 --beams pp:10TeV -P pythia6.params --randomize-seed -n $nevt -o ./events.fifo &
time rivet -a FASTSIM -H fastsim -x 1 -n ${nevt}  -l Rivet.Analysis=${LEVEL} events.fifo
#hepmc


