# Folder for generated mod-files
set(moduledir "${CMAKE_CURRENT_BINARY_DIR}/modules")

add_library(scalapackfx)
add_library(scalapackfx::scalapackfx ALIAS scalapackfx)

set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

set(sources-fpp
    ${curdir}/blacs.fpp
    ${curdir}/blacsfx.fpp
    ${curdir}/blacsgrid.fpp
    ${curdir}/libscalapackfx.fpp
    ${curdir}/linecomm.fpp
    ${curdir}/pblas.fpp
    ${curdir}/pblasfx.fpp
    ${curdir}/scalapack.fpp
    ${curdir}/scalapackfx.fpp
    ${curdir}/scalapackfx_common.fpp
    ${curdir}/scalapackfx_tools.fpp
)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp})

scalapackfx_preprocess("${FYPP}" "${fypp_flags}" "fpp" "f90" "${ALL-SOURCES-FPP}" all-sources-f90-preproc)

target_sources(
  scalapackfx
  PRIVATE
  ${all-sources-f90-preproc}
)

target_include_directories(
  scalapackfx PUBLIC
  $<BUILD_INTERFACE:${moduledir}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/${SCALAPACKFX_INSTALL_MODULEDIR}>
)

target_link_libraries(scalapackfx PRIVATE MPI::MPI_Fortran)

set_target_properties(
  scalapackfx
  PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    EXPORT_NAME scalapackfx
    OUTPUT_NAME scalapackfx
    Fortran_MODULE_DIRECTORY "${moduledir}"
)

if (SCALAPACKFX_INSTALL)
  install(
    TARGETS scalapackfx
    EXPORT scalapackfxTargets
    LIBRARY
      DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      COMPONENT scalapackfx_runtime
      NAMELINK_COMPONENT scalapackfx_development
    PUBLIC_HEADER
      DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
      COMPONENT scalapackfx_development
  )
  install(
    DIRECTORY "${moduledir}/"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/${SCALAPACKFX_INSTALL_MODULEDIR}"
    COMPONENT scalapackfx_development
  )
endif ()
