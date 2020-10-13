set(sources-fpp
  blacs.fpp
  blacsfx.fpp
  blacsgrid.fpp
  libscalapackfx.fpp
  linecomm.fpp
  pblas.fpp
  pblasfx.fpp
  scalapack.fpp
  scalapackfx.fpp
  scalapackfx_common.fpp
  scalapackfx_tools.fpp)

fypp_preprocess("${sources-fpp}" sources-f90)

add_library(scalapackfx ${sources-f90})

target_link_libraries(scalapackfx PRIVATE Scalapack::Scalapack LAPACK::LAPACK)

set(BUILD_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/include)

set_target_properties(scalapackfx PROPERTIES Fortran_MODULE_DIRECTORY ${BUILD_MOD_DIR})

target_include_directories(scalapackfx PUBLIC
  $<BUILD_INTERFACE:${BUILD_MOD_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${INSTALL_MODULEDIR}>)

install(TARGETS scalapackfx
  EXPORT scalapackfx-targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

if(INSTALL_INCLUDE_FILES)
  install(DIRECTORY ${BUILD_MOD_DIR}/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${INSTALL_MODULEDIR})
endif()
