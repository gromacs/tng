if(NOT OPENMP_FOUND)
  add_library(openmp SHARED openmp_stubs.c)
  set(OpenMP_LIBS "openmp")
endif()

add_library(tng_io SHARED tng_io.c md5.c)

