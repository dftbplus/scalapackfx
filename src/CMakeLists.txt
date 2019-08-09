set(sources-fpp
  blacs.F90
  blacsfx.F90
  blacsgrid.F90
  libscalapackfx.F90
  linecomm.F90
  pblas.F90
  pblasfx.F90
  scalapack.F90
  scalapackfx.F90
  scalapackfx_common.F90
  scalapackfx_tools.F90)

foreach(fppsrc IN LISTS sources-fpp)
  string(REGEX REPLACE "\\.F90" ".f90" f90src ${fppsrc})  
  string(REGEX REPLACE "\\.F90" ".m4" m4src ${fppsrc})
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${m4src})
    set(otherdeps ${CMAKE_CURRENT_SOURCE_DIR}/${m4src})
  else()
    set(otherdeps)
  endif()
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${f90src}
    COMMAND m4 -I${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${fppsrc} > ${CMAKE_CURRENT_BINARY_DIR}/${f90src}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${fppsrc} ${otherdeps})
  list(APPEND SOURCES-F90-PREPROC ${CMAKE_CURRENT_BINARY_DIR}/${f90src})
endforeach()

add_library(scalapackfx_objlib OBJECT ${SOURCES-F90-PREPROC})
add_library(scalapackfx $<TARGET_OBJECTS:scalapackfx_objlib>)

set(includedir ${CMAKE_CURRENT_BINARY_DIR}/include)

set_target_properties(scalapackfx_objlib PROPERTIES Fortran_MODULE_DIRECTORY ${includedir})

target_include_directories(scalapackfx_objlib PUBLIC ${includedir})
target_include_directories(scalapackfx PUBLIC ${includedir})

install(TARGETS scalapackfx
  ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
  LIBRARY DESTINATION ${INSTALL_LIB_DIR})

install(DIRECTORY ${includedir}/ DESTINATION ${INSTALL_MOD_DIR})
