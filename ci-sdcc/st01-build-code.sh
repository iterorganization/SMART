#!/bin/bash
# Bamboo CI script to build code and create library
# Execute script from root directory

# setup environment
# Get toolchain version
source ./ci-sdcc/st00-header.sh $1 $2

# Note Disable set -e option when using on local as it will exit the shell on error

if [[ "$(uname -n)" == *"bamboo"* ]]; then
   set -e -u -o pipefail
fi
LIBRARY_NAME=libsmart

echo "Compiling code"
set -x
make clean validate lib exe FC="$FCOMPILER"
set +x
echo "Done compiling code"
libfilepath="./$LIBRARY_NAME.a"
includefilepath="./smart.mod"

if [ ! -f "$libfilepath" ]; then
    echo "$LIBRARY_NAME.a does not exist: $libfilepath"
    # exit 1
fi

ARTIFACT="lib.tar.gz"
echo "Checking if artifact exists..."
if [ -f "$ARTIFACT" ]; then
    rm "$ARTIFACT"
    echo "$ARTIFACT removed successfully."
fi

# Create acrtifact
echo "Creating artifact..."
tar -cvzf $ARTIFACT "$libfilepath" "$includefilepath" >/dev/null 2>&1
if [ -f "$ARTIFACT" ]; then
    echo "Artifact $ARTIFACT created successfully."
fi

# show contents of artifact
echo "Showing contents of Artifact"
tar -tzvf $ARTIFACT

# Cleanup

echo "Done"
