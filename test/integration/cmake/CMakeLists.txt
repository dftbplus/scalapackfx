cmake_minimum_required(VERSION 3.16)

project(TestScalapackFxBuild LANGUAGES Fortran)

find_package(ScalapackFx REQUIRED)

add_executable(test_build test_build.f90)
target_link_libraries(test_build ScalapackFx::ScalapackFx)
