#!/bin/bash
BASEDIR=$HOME/local
#HEPMCDIR=/afs/cern.ch/sw/lcg/external/HepMC/2.04.01/slc4_amd64_gcc34/include
#BOOSTDIR=/afs/cern.ch/sw/lcg/external/Boost/1.37.0_python2.5/slc4_amd64_gcc34/include/boost-1_37/

#ATHENASW=/home/disipio/datanas/athena/15.1.0/sw/lcg/external

HEPMCDIR=$BASEDIR/include
#BOOSTDIR=$ATHENASW/Boost/1.38.0_python2.5/slc4_ia32_gcc34/include/boost-1_38
BOOSTDIR=/home/disipio/local/boost/include/boost-1_37

     #-I$ATHENASW/MCGenerators/rivet/1.1.1/slc4_ia32_gcc34/include \


#CFLAGS='-m32'
CFLAGS=''

if [ "`uname`" == "Darwin" ]
then
 LIBEXT=dylib
 SHAREDOPT='-dynamiclib -flat_namespace -undefined dynamic_lookup -single_module'
else 
 LIBEXT=so
 SHAREDOPT='-shared -fPIC -DPIC'
fi
echo Compiling for platform `uname`

target="fastsim"
while [ $# -gt 0 ] ; do
  case $1 in
  "OPT")
    SYMBOL="-O2 $SYMBOL"
    shift 1
    ;;
  "ALLOPT")
    SYMBOL="-D FWM -D YLM $SYMBOL"
    shift 1
    ;;
  "DEBUG")
    SYMBOL="-O0 -g $SYMBOL"
    shift 1
    ;;
  *)
    target="fastsim"
    shift 1
    ;;
  esac
done

if [[ ! -z $SYMBOL ]] ; then
  echo "Compiling with defined symbol $SYMBOL"
fi

echo "Compiling $target"
case $target in
  "fastsim")
     g++ -Wall ${CFLAGS} -o libRivetDetectorFastSim.$LIBEXT $SHAREDOPT $SYMBOL -I$BASEDIR/include \
         -I$HEPMCDIR -I$BOOSTDIR -I`root-config --incdir` `root-config --libs` `root-config --glibs` DetectorFastSim.cc
     ;;
    *)
     echo "Unknown target $target"
     exit 1
     ;;
esac
    
if [[ $? -ne 0 ]]
then
  echo Ok
fi
