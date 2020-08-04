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

foreach(fppsrc IN LISTS sources-fpp)
  string(REGEX REPLACE "\\.fpp" ".f90" f90src ${fppsrc})  
  string(REGEX REPLACE "\\.fpp" ".fypp" fyppsrc ${fppsrc})
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${fyppsrc})
    set(otherdeps ${CMAKE_CURRENT_SOURCE_DIR}/${fyppsrc})
  else()
    set(otherdeps)
  endif()
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${f90src}
    COMMAND ${FYPP} -I${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${fppsrc} ${CMAKE_CURRENT_BINARY_DIR}/${f90src}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${fppsrc} ${otherdeps})
  list(APPEND SOURCES-F90-PREPROC ${CMAKE_CURRENT_BINARY_DIR}/${f90src})
endforeach()

add_library(scalapackfx ${SOURCES-F90-PREPROC})
target_link_libraries(scalapackfx ${EXTERNAL_LIBRARIES} ${MPI_Fortran_LIBRARIES})

set(BUILD_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/include)

set_target_properties(scalapackfx PROPERTIES Fortran_MODULE_DIRECTORY ${BUILD_MOD_DIR})

target_include_directories(scalapackfx PUBLIC
  $<BUILD_INTERFACE:${BUILD_MOD_DIR}>
  $<INSTALL_INTERFACE:${INSTALL_MOD_DIR}>)

install(TARGETS scalapackfx
  EXPORT ${INSTALL_EXPORT_NAME}
  ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
  LIBRARY DESTINATION ${INSTALL_LIB_DIR})

if(INSTALL_INCLUDE_FILES)
  install(DIRECTORY ${BUILD_MOD_DIR}/ DESTINATION ${INSTALL_MOD_DIR})
endif()
