#!/bin/bash
#
# Tests whether the installed ScalapackFx library can be used within a CMake project.
#
# Arguments:
#
#   - building directory (will be created, should not exist)
#
# Requirements:
#
#   - Environment variable FC contains the same Fortran compiler as used for ScalapackFx
#
#   - Environment variable CMAKE_PREFIX_PATH contains the ScalapackFx install root.
#
SCRIPTDIR=$(dirname $0)
SCRIPTNAME=$(basename $0)
BUILDDIR=$1
shift

if [ -d ${BUILDDIR} ]; then
  echo "${SCRIPTNAME}: Test build directory '${BUILDDIR}' already exists." >&2
  exit 1
fi

FC=$FC cmake -B ${BUILDDIR} ${SCRIPTDIR} -DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH} "$@" \
  || { echo "Configuration step failed" >&2; exit 1; }
cmake --build ${BUILDDIR} -- VERBOSE=1 || { echo "Build step failed" >&2; exit 1; }
echo "CMake build succeeded!"
