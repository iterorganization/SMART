#!/bin/bash
source /etc/profile.d/modules.sh
source ./ci-sdcc/utils.sh
##########################################################################################
#                     Set environment based on toolchain                                 #
##########################################################################################
module use /work/imas/etc/modules/all
module use /work/imas/etc/modules-2025b/all 2>/dev/null || true

# expand aliases
shopt -s expand_aliases

#print hostname
hostname -f

# Pick up toolchain from $1 (default foss-2023b).  $2 is preserved for
# backward compatibility with existing Bamboo plan calls but is no longer
# consulted: the dep files below pin a specific AL/HLI version per cell.
if [ -n "$1" ]; then
    echo "> Compiling with $1.  Previously loaded modules will be purged."
    module purge
    TOOLCHAIN_VERSION="$1"
fi

if [ -z "$TOOLCHAIN_VERSION" ]; then
    echo "> No toolchain found, setting it to default : foss-2023b"
    TOOLCHAIN_VERSION="foss-2023b"
fi

echo "> Building for $TOOLCHAIN_VERSION"

if [[ $TOOLCHAIN_VERSION == *"intel"* ]]; then
    FC="ifort"
fi
if [[ $TOOLCHAIN_VERSION == *"foss"* ]]; then
    FC="gfortran"
fi

buildtime_dependencies="./ci-sdcc/deps/${TOOLCHAIN_VERSION}/buildtime.txt"
runtime_dependencies="./ci-sdcc/deps/${TOOLCHAIN_VERSION}/runtime.txt"
if [ ! -f "$buildtime_dependencies" ]; then
    echo "File $buildtime_dependencies not found."
    return 1
fi
if [ ! -f "$runtime_dependencies" ]; then
    echo "File $runtime_dependencies not found."
    return 1
fi
echo "> Listing available modules"
echo "-------------------------------------------------------"
echo "> build time modules"

declare -a BUILDMODULES=()
declare -a RUNMODULES=()
declare -a EBBUILDMODULES=()
declare -a EBBRUNMODULES=()
counter=0
# Read the buildtime file.  Each non-empty line must be a fully-qualified
# module name (containing a '/'); bare basenames fall back to getModuleName
# in utils.sh, but are discouraged.
while IFS= read -r line || [[ -n $line ]]; do
    if [[ -z "${line// /}" ]]; then
        counter=$(("$counter" + 1))
        continue
    fi
    if [[ $line == *"/"* ]]; then
        echo "  fixed module: $line"
        BUILDMODULES["$counter"]="$line"
        EBBUILDMODULES["$counter"]="('$line', EXTERNAL_MODULE),"
    else
        module_version=$(getModuleName "$line" "$TOOLCHAIN_VERSION" "")
        echo "  resolved $line -> $module_version"
        BUILDMODULES["$counter"]="$module_version"
        EBBUILDMODULES["$counter"]=$(getModuleNameAndVersion "$module_version")
    fi
    counter=$(("$counter" + 1))
done <"$buildtime_dependencies"
echo "-------------------------------------------------------"
echo "> run time modules"

counter=0
while IFS= read -r line || [[ -n $line ]]; do
    line="${line// /}"
    if [[ -z "$line" ]]; then
        counter=$(("$counter" + 1))
        continue
    fi
    if [[ $line == *"/"* ]]; then
        echo "  fixed module: $line"
        RUNMODULES["$counter"]="$line"
        EBBRUNMODULES["$counter"]="('$line', EXTERNAL_MODULE),"
    else
        module_version=$(getModuleName "$line" "$TOOLCHAIN_VERSION" "")
        echo "  resolved $line -> $module_version"
        RUNMODULES["$counter"]="$module_version"
        EBBRUNMODULES["$counter"]=$(getModuleNameAndVersion "$module_version")
    fi
    counter=$(("$counter" + 1))
done <"$runtime_dependencies"
echo "-------------------------------------------------------"

echo "> Details of environment"
echo "    TOOLCHAIN_VERSION : $TOOLCHAIN_VERSION"
echo "    BUILDMODULES :   ${BUILDMODULES[*]}"
echo "    RUNMODULES :     ${RUNMODULES[*]}"
echo "    EBBUILDMODULES : ${EBBUILDMODULES[*]}"
echo "    EBRUNMODULES :   ${EBBRUNMODULES[*]}"
echo "    Compiler : $FC"
echo "-------------------------------------------------------"

echo "> Loading build time modules"
for imodule in "${BUILDMODULES[@]}"; do
    IFS='/' read -r iname iversion <<< "$imodule"
    MODULE_EXISTS=$(module -r -t list 2>&1 | grep -E "$iname/")
    if [ -z "$MODULE_EXISTS" ]; then
        echo "    $iname not available, loading $imodule"
        module load "$imodule"
    else
        echo "    $MODULE_EXISTS already loaded"
    fi
done

echo "> Loading run time modules"
for imodule in "${RUNMODULES[@]}"; do
    IFS='/' read -r iname iversion <<< "$imodule"
    MODULE_EXISTS=$(module -r -t list 2>&1 | grep -E "$iname/")
    if [ -z "$MODULE_EXISTS" ]; then
        echo "    $iname not available, loading $imodule"
        module load "$imodule"
    else
        echo "    $MODULE_EXISTS already loaded"
    fi
done

echo "> Done"
echo "-------------------------------------------------------"
