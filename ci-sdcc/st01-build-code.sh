#!/usr/bin/env bash
# Bamboo CI script to build code and create library
# Execute script from root directory

# Capture positional args before sourcing st00 (which sets $1/$2 during
# its execution).  $1 = toolchain (foss-2023b, intel-2023b, ...), $2 = AL
# major version (5).  The artifact filename below embeds both because the
# Bamboo plan publishes them as e.g. lib-foss-2023b-al5.
TOOLCHAIN="${1:-foss-2023b}"
AL_MAJOR="${2:-5}"

# setup environment
# Get toolchain version
source ./ci-sdcc/st00-header.sh "$TOOLCHAIN" "$AL_MAJOR"

# Note Disable set -e option when using on local as it will exit the shell on error

if [[ "$(uname -n)" == *"bamboo"* ]]; then
   set -e -u -o pipefail
fi
LIBRARY_NAME=libsmart

echo "Compiling code"
set -x
make clean validate lib exe FC="$FC"
set +x
echo "Done compiling code"
libfilepath="./$LIBRARY_NAME.a"
includefilepath="./mod_smart.mod"

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
