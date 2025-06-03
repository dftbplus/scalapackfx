# Variables influencing how subproject is obtained
set(CMAKE_REQUIRE_FIND_PACKAGE_Fortuno ${SCALAPACKFX_SUBPROJECT_REQUIRE_FIND})
set(CMAKE_DISABLE_FIND_PACKAGE_Fortuno ${SCALAPACKFX_SUBPROJECT_DISABLE_FIND})
# set FETCHCONTENT_SOURCE_DIR_FORTUNO to use a local source of the subproject

# Subproject related variables
option(
  FORTUNO_BUILD_SHARED_LIBS "Fortuno: Build as shared library" ${SCALAPACKFX_BUILD_SHARED_LIBS}
)

option(FORTUNO_WITH_MPI "Fortuno: whether to build the MPI interface" ON)

# Make subproject available
FetchContent_Declare(
  Fortuno
  GIT_REPOSITORY "https://github.com/fortuno-repos/fortuno.git"
  GIT_TAG "main"
  FIND_PACKAGE_ARGS
)
FetchContent_MakeAvailable(Fortuno)

if (Fortuno_FOUND)
  message(STATUS "Subproject Fortuno: using installed version")
else ()
  message(STATUS "Subproject Fortuno: building from source in ${fortuno_SOURCE_DIR}")
endif ()
