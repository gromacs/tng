if(NOT OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

include_directories(${TRAJECTORY_SOURCE_DIR}/src/lib)
link_directories(${TRAJECTORY_BINARY_DIR}/src/lib)

set(EXAMPLE_FILES_DIR ${CMAKE_BINARY_DIR}/example_files/ CACHE STRING "Directory where to find input test files and save output files.")

configure_file(tng_io_testing.h.in ${CMAKE_BINARY_DIR}/generated/tng_io_testing.h)
include_directories(${CMAKE_BINARY_DIR}/generated/)

add_executable(tng_testing tng_io_testing.c)
target_link_libraries(tng_testing tng_io)

add_executable(md_openmp md_openmp.c)
target_link_libraries(md_openmp tng_io ${OpenMP_LIBS} m)