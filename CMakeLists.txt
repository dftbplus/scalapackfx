cmake_minimum_required(VERSION 3.16)

include(CMakePackageConfigHelpers)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)
include(ScalapackFxUtils)

include(${CMAKE_CURRENT_SOURCE_DIR}/config.cmake)

project(ScalapackFx VERSION 1.0.2 LANGUAGES Fortran)

setup_build_type()

#
# Prerequisites
#
find_package(CustomScalapack REQUIRED QUIET)
find_package(CustomLapack REQUIRED)
find_program(FYPP fypp)
if(NOT FYPP)
  message(FATAL_ERROR "Preprocessor fypp could not be found")
endif()

#
# Build instructions
#
include(GNUInstallDirs)

add_subdirectory(lib)
if(NOT BUILD_EXPORTED_TARGETS_ONLY)
  add_subdirectory(test)
endif()

#
# Installation
#
add_library(ScalapackFx INTERFACE)
target_link_libraries(ScalapackFx INTERFACE scalapackfx)
install(TARGETS ScalapackFx EXPORT scalapackfx-targets)

install(EXPORT scalapackfx-targets
  FILE scalapackfx-targets.cmake
  NAMESPACE ScalapackFx::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/scalapackfx"
  EXPORT_LINK_INTERFACE_LIBRARIES)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/utils/export/scalapackfx-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/scalapackfx-config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scalapackfx)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/scalapackfx-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/scalapackfx-config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/scalapackfx-config-version.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scalapackfx)

install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scalapackfx)

get_pkgconfig_params(PKGCONFIG_REQUIRES PKGCONFIG_LIBS PKGCONFIG_LIBS_PRIVATE PKGCONFIG_C_FLAGS)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/utils/export/scalapackfx.pc.in
  ${CMAKE_CURRENT_BINARY_DIR}/scalapackfx.pc @ONLY)
install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/scalapackfx.pc"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
