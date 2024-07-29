#!/bin/bash
# Bamboo CI script to build actor and run standalone program
# Execute script from root directory

source ./ci-sdcc/st00-header.sh $1 $2
USERNAME=$(whoami)
# Note Disable set -e option when using on local as it will exit the shell on error
if [[ "$(uname -n)" == *"bamboo"* ]]; then
    set -e -u -o pipefail
fi
ACTOR_NAME=smart
export ACTOR_FOLDER=./$ACTOR_NAME
# remove actor directory if present
if [ -d $ACTOR_FOLDER ]; then
    rm -rf $ACTOR_FOLDER
    echo "$ACTOR_FOLDER removed successfully."
fi

set -x
# create actor
make actor FC="$FCOMPILER"

# Run actor standalone program
# python run_hcd2core_profilescode --src "imas:mdsplus?user=public;pulse=130012;run=115;database=TEST;version=3" --dest "imas:mdsplus?user=$USERNAME;pulse=130012;run=23;database=ITER;version=3" --time 200.0
set +x

# remove __pycache__ from directory
find $ACTOR_FOLDER -type d -name '__pycache__' -exec rm -rf {} +

# Create acrtifact
tar -cvzf actor.tar.gz $ACTOR_NAME >/dev/null 2>&1

set -x
# show contents of artifact
tar -tzvf actor.tar.gz

set +x

if [[ "$(uname -n)" != "sdcc"* ]]; then
    rm -rf $ACTOR_FOLDER
fi

echo "Done"
