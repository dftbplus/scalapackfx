# Folder for generated mod-files
set(moduledir "${CMAKE_CURRENT_BINARY_DIR}/modules")

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
  set_target_properties(
    ${test} PROPERTIES Fortran_MODULE_DIRECTORY "${moduledir}"
  )
  target_link_libraries(${test} PRIVATE scalapackfx)
  target_link_libraries(${test} PRIVATE Fortuno::fortuno_mpi MPI::MPI_Fortran)
  target_link_libraries(${test} PRIVATE Scalapack::Scalapack)

  add_test(
    NAME ${test}
    COMMAND ${test}
  )
endforeach()
