# Folder for generated mod-files
set(moduledir "${CMAKE_CURRENT_BINARY_DIR}/modules")

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
#  test_tran)

add_executable(test_desc)
target_sources(test_desc PRIVATE test_desc.f90)
set_target_properties(
  test_desc PROPERTIES Fortran_MODULE_DIRECTORY "${moduledir}"
)
target_link_libraries(test_desc PRIVATE scalapackfx)
target_link_libraries(test_desc PRIVATE Fortuno::fortuno_mpi MPI::MPI_Fortran)
target_link_libraries(test_desc PRIVATE Scalapack::Scalapack)

add_test(
  NAME test_desc
  COMMAND test_desc
)
