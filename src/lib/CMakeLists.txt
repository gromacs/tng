if(TNG_BUILD_FORTRAN)
  add_library(tng_io tng_io.c md5.c tng_io_fortran.c)
else()
  add_library(tng_io tng_io.c md5.c)
endif()

set_property(TARGET tng_io PROPERTY LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

install(TARGETS tng_io
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

if(HAVE_INTTYPES_H)
  set_property(TARGET tng_io APPEND PROPERTY COMPILE_DEFINITIONS USE_STD_INTTYPES_H)
endif()

if(ZLIB_FOUND)
  set_property(TARGET tng_io APPEND PROPERTY COMPILE_DEFINITIONS USE_ZLIB)
  target_link_libraries(tng_io tng_compress ${ZLIB_LIBRARIES})
else()
  target_link_libraries(tng_io tng_compress)
endif()
