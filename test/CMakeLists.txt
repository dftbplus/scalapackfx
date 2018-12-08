set(TARGETS
  test_desc
  test_subgrids)

set(COMMON-DEP-TARGETS
  test_det
  test_diag
  test_psyr_pher
  test_svd
  test_tran)

foreach(target ${TARGETS})
  add_executable(${target} ${target}.f90)
  target_link_libraries(${target} scalapackfx ${SCALAPACK_LIBRARIES} ${MPI_Fortran_LIBRARIES})
endforeach()

foreach(target ${COMMON-DEP-TARGETS})
  add_executable(${target} ${target}.f90 test_common.f90)
  target_link_libraries(${target} scalapackfx ${SCALAPACK_LIBRARIES} ${MPI_Fortran_LIBRARIES})
endforeach()
