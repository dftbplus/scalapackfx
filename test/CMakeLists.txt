# Folder for generated mod-files
set(moduledir "${CMAKE_CURRENT_BINARY_DIR}/modules")

add_library(blacstestutils blacstestutils.f90)
target_link_libraries(blacstestutils PRIVATE scalapackfx)
target_link_libraries(blacstestutils PRIVATE Fortuno::fortuno_mpi)

list(APPEND tests
  test_desc
#  test_cpg2l
#  test_det
#  test_diag
#  test_gemr2d
#  test_linecomm
#  test_matinv
#  test_pposv
#  test_psyr_pher
#  test_remoteelements
#  test_svd
#  test_tran
)


foreach(test IN LISTS tests)
  add_executable(${test})
  target_sources(${test} PRIVATE ${test}.f90)
  #set_target_properties(
  #  ${test} PROPERTIES Fortran_MODULE_DIRECTORY "${moduledir}"
  #)
  target_link_libraries(${test} PRIVATE blacstestutils)
  target_link_libraries(${test} PRIVATE scalapackfx)
  target_link_libraries(${test} PRIVATE Fortuno::fortuno_mpi MPI::MPI_Fortran)
  target_link_libraries(${test} PRIVATE Scalapack::Scalapack)
  add_test(
    NAME ${test}
    COMMAND  ${MPI_TEST_RUNNER} ./${test}
  )
  set_tests_properties(
    ${test}
    PROPERTIES ENVIRONMENT OMP_NUM_THREADS=${SCALAPACKFX_TEST_OMP_THREADS}
  )
endforeach()
