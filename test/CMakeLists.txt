set(targets
  test_desc
  test_subgrids)

set(common-dep-targets
  test_det
  test_diag
  test_linecomm
  test_psyr_pher
  test_svd
  test_tran)

add_library(common OBJECT test_common.f90)
#target_link_libraries(common PRIVATE scalapackfx)
add_dependencies(common scalapackfx)
target_include_directories(common
  PRIVATE $<TARGET_PROPERTY:scalapackfx,INTERFACE_INCLUDE_DIRECTORIES>)

foreach(target IN LISTS targets)
  add_executable(${target} ${target}.f90)
  target_link_libraries(${target} scalapackfx ${EXTERNAL_LIBRARIES} ${MPI_Fortran_LIBRARIES})
endforeach()

foreach(target IN LISTS common-dep-targets)
  add_executable(${target} ${target}.f90 $<TARGET_OBJECTS:common>)
  target_link_libraries(${target} scalapackfx ${EXTERNAL_LIBRARIES}
    ${MPI_Fortran_LIBRARIES})
endforeach()
