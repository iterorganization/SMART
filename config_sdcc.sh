# Start from clean environment
ml purge

# IMAS and iWrap
# intel
ml IMAS-AL-Fortran/5.2.2-intel-2023b-DD-3.41.0
ml IMAS-AL-Python/5.2.2-intel-2023b-DD-3.41.0
ml XMLlib/3.3.2-intel-compilers-2023.2.1
# foss
#ml IMAS-AL-Fortran/5.2.2-foss-2023b-DD-3.41.0
#ml IMAS-AL-Python/5.2.2-foss-2023b-DD-3.41.0
#ml XMLlib/3.3.2-GCC-13.2.0

# iWrap
ml iWrap/0.10.0-GCCcore-13.2.0

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
