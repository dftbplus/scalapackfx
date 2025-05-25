#[=============================================================================[
#                       Basic project definition                              #
]=============================================================================]

cmake_minimum_required(VERSION 3.22...3.29)

list(APPEND CMAKE_MESSAGE_CONTEXT scalapackfx)

project(
  scalapackfx
  VERSION 1.3.0
  DESCRIPTION "Modern Fortran Interface for ScaLAPACK"
  LANGUAGES Fortran
)

#[=============================================================================[
#                       Options                                               #
]=============================================================================]

option(
  SCALAPACKFX_BUILD_SHARED_LIBS
  "scalapackfx: Build as shared library"
  OFF
)

option(
  SCALAPACKFX_BUILD_EXAMPLES
  "scalapackfx: Whether to build the examples"
  ${PROJECT_IS_TOP_LEVEL}
)
option(
  SCALAPACKFX_BUILD_TESTS
  "scalapackfx: Whether to build the tests"
  ${PROJECT_IS_TOP_LEVEL}
)
option(
  SCALAPACKFX_INSTALL
  "scalapackfx: Install project"
  ${PROJECT_IS_TOP_LEVEL}
)
set(
  SCALAPACKFX_INSTALL_MODULEDIR
  "modules"
  CACHE STRING
  "scalapackfx: Sub-directory for installed Fortran module files (relative to CMAKE_INSTALL_LIBDIR)"
)

option(
  SCALAPACKFX_SUBPROJECT_REQUIRE_FIND
  "scalapackfx: Require find_package for all subprojects"
  OFF
)
option(
  SCALAPACKFX_SUBPROJECT_DISABLE_FIND
  "scalapackfx: Disable find_package for all subprojects"
  OFF
)

#[=============================================================================[
#                  Project configuration                                      #
]=============================================================================]

include(FetchContent)
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)
include(scalapackfxHelpers)
include(ScalapackFxUtils)

set(BUILD_SHARED_LIBS ${SCALAPACKFX_BUILD_SHARED_LIBS})
scalapackfx_setup_build_type("RelWithDebInfo")

find_package(MPI REQUIRED)

#
# Prerequisites
#
find_package(CustomScalapack REQUIRED QUIET)
find_package(CustomLapack REQUIRED)

include(subprojects/fypp.cmake)

# Extra flags from cmake options and compiler capabilities
scalapackfx_add_fypp_defines(FYPP_FLAGS)

set(FYPP_FLAGS "${FYPP_FLAGS}")

set(FYPP_CONFIG_FLAGS "${FYPP_FLAGS}")
# Make sure, the line-marker option is not set
list(REMOVE_ITEM FYPP_CONFIG_FLAGS "-n")
set(FYPP_BUILD_FLAGS "${FYPP_FLAGS}" "--file-var-root=${CMAKE_SOURCE_DIR}"
  "$<IF:$<CONFIG:Debug>,-DDEBUG=1,-DDEBUG=0>")

set(PYTHON_INTERPRETER "python3" CACHE STRING
  "Python interpreter to use for preprocessor")

set(fypp_flags ${FYPP_BUILD_FLAGS})
list(APPEND
  fypp_flags -I${CMAKE_CURRENT_SOURCE_DIR}/include -DRELEASE="'${VERSION}'"
)

#[=============================================================================[
#                       Main definition                                       #
]=============================================================================]

add_subdirectory(src)

if (SCALAPACKFX_BUILD_EXAMPLES)
  add_subdirectory(example)
endif ()
if (SCALAPACKFX_BUILD_TESTS)
  include(subprojects/Fortuno.cmake)
  enable_testing()
  add_subdirectory(test)
endif ()


#[=============================================================================[
#                       Install or Export                                     #
]=============================================================================]

if (SCALAPACKFX_INSTALL)

  # pkg-config files
  configure_file(cmake/scalapackfx.pc.in scalapackfx.pc @ONLY)
  install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/scalapackfx.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
    COMPONENT scalapackfx_development
  )

  # cmake export files
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/scalapackfxConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
  )
  configure_package_config_file(
      cmake/scalapackfxConfig.cmake.in
      ${CMAKE_CURRENT_BINARY_DIR}/scalapackfxConfig.cmake
      INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scalapackfx
  )
  install(
    FILES
      ${CMAKE_CURRENT_BINARY_DIR}/scalapackfxConfigVersion.cmake
      ${CMAKE_CURRENT_BINARY_DIR}/scalapackfxConfig.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scalapackfx
    COMPONENT scalapackfx_development
  )

  export(
    EXPORT scalapackfxTargets
    FILE scalapackfxTargets.cmake
    NAMESPACE scalapackfx::
  )
  install(
    EXPORT scalapackfxTargets
    FILE scalapackfxTargets.cmake
    NAMESPACE scalapackfx::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scalapackfx
    COMPONENT scalapackfx_development
  )
endif ()

# Make project available for FetchContent
if (NOT PROJECT_IS_TOP_LEVEL)
  # Propagate variables
  set(scalapackfx_VERSION ${scalapackfx_VERSION} PARENT_SCOPE)
  set(scalapackfx_VERSION_MAJOR ${scalapackfx_VERSION_MAJOR} PARENT_SCOPE)
  set(scalapackfx_VERSION_MINOR ${scalapackfx_VERSION_MINOR} PARENT_SCOPE)
  set(scalapackfx_VERSION_PATCH ${scalapackfx_VERSION_PATCH} PARENT_SCOPE)
  set(scalapackfx_VERSION_TWEAK ${scalapackfx_VERSION_TWEAK} PARENT_SCOPE)
endif ()
