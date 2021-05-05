cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
enable_language(Fortran)
project(TALYS VERSION 1.95 LANGUAGES Fortran)


# First step is to configure the machine file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/source/machine.f.in ${CMAKE_CURRENT_SOURCE_DIR}/source/machine.f)

# Now we can find all the source files
file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/source/*.f")


# Last we can make the binary
add_executable(talys ${sources})

