cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
enable_language(Fortran)
project(TALYS VERSION 2.0 LANGUAGES Fortran)

include(GNUInstallDirs)


# First step is to configure the machine file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/source/machine.f90.in ${CMAKE_CURRENT_SOURCE_DIR}/source/machine.f90)

# Now we can find all the source files
file(GLOB_RECURSE sources_f CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/source/*.f")
file(GLOB_RECURSE sources_f90 CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/source/*.f90")

# Last we can make the binary
add_executable(talys ${sources_f} ${sources_f90})

INSTALL(TARGETS talys RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})