add_library(tng_io tng_io.c md5.c)

if(HAVE_INTTYPES_H)
  set_property(TARGET tng_io APPEND PROPERTY COMPILE_DEFINITIONS USE_STD_INTTYPES_H)
endif()

if(BUILD_FORTRAN)
  set_property(TARGET tng_io APPEND PROPERTY COMPILE_DEFINITIONS BUILD_FORTRAN)
endif()

if(ZLIB_FOUND)
  set_property(TARGET tng_io APPEND PROPERTY COMPILE_DEFINITIONS USE_ZLIB)
  target_link_libraries(tng_io tng_compress ${ZLIB_LIBRARIES})
else()
  target_link_libraries(tng_io tng_compress)
endif()
