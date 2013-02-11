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

if(OPENMP_FOUND)
  add_executable(md_openmp md_openmp.c)
  target_link_libraries(md_openmp tng_io ${OpenMP_LIBS} m)
endif()

add_executable(tng_io_read_pos tng_io_read_pos.c)
target_link_libraries(tng_io_read_pos tng_io)

add_executable(tng_parallel_read tng_parallel_read.c)
target_link_libraries(tng_parallel_read tng_io)

if(BUILD_FORTRAN)
  # This does not work due to a bug in CMake. Remove lines below if no fortran compiler is found.
  enable_language(Fortran OPTIONAL)
  if(${CMAKE_Fortran_COMPILER_WORKS})
    get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
    if (Fortran_COMPILER_NAME STREQUAL "gfortran")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcray-pointer ${OpenMP_C_FLAGS} -std=legacy")
    endif()
    if(OPENMP_FOUND)
      add_executable(md_openmp_f md_openmp.f)
      target_link_libraries(md_openmp_f tng_io ${OpenMP_LIBS})
    endif()
  endif()
endif()
