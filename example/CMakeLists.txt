# Folder for generated mod-files
set(moduledir "${CMAKE_CURRENT_BINARY_DIR}/modules")

add_executable(transpose)

target_sources(transpose PRIVATE transpose.f90 common_module.f90)
set_target_properties(transpose
  PROPERTIES Fortran_MODULE_DIRECTORY "${moduledir}"
)
target_link_libraries(transpose PRIVATE scalapackfx)
target_link_libraries(transpose PRIVATE MPI::MPI_Fortran)
target_link_libraries(transpose PRIVATE Scalapack::Scalapack)
