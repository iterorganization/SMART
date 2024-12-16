# Start from clean environment

#COMPILER:ifx/gfortran

export COMPILER=ifx
#export COMPILER=gfortran

# IMAS and iWrap
# intel
if [[ $COMPILER =~ ^(ifort|icc|icpc|intel|ifx|icx)$ ]]
then
  ml purge
  ml IMAS-AL-Fortran/5.3.0-intel-2023b-DD-3.42.0
  ml IMAS-AL-Python/5.3.0-intel-2023b-DD-3.42.0
  ml XMLlib/3.3.2-intel-compilers-2023.2.1
# foss
elif [[ $COMPILER =~ ^(gfortran|g++|gcc|GCC)$ ]]
then
  ml purge
  ml IMAS-AL-Fortran/5.3.0-foss-2023b-DD-3.42.0
  ml IMAS-AL-Python/5.3.0-foss-2023b-DD-3.42.0
  ml XMLlib/3.3.2-GCC-13.2.0
#
else
  echo "Set environt variable as ifx/gfortran"
  return
fi

# Constants
ml Fundamental-Constants

# iWrap
ml iWrap/1.0.0-GCCcore-13.2.0

# For waveform editions
#module load Waveform-Cooker

# For debugging, just in case
#ml TotalView


# Actor folder
export ACTOR_FOLDER=~/public/PYTHON_ACTORS
mkdir -p $ACTOR_FOLDER

# Need to remove the stack limit to avoid segmentation fault inside codes
ulimit -Ss unlimited

# EXTEND PYTHON PATH AND AVOID DOUBLONS
export PYTHONPATH=$ACTOR_FOLDER:$PYTHONPATH
export PYTHONPATH="$(perl -e 'print join(":", grep { not $seen{$_}++ } split(/:/, $ENV{PYTHONPATH}))')"
