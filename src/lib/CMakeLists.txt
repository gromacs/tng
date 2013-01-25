add_library(tng_io SHARED tng_io.c md5.c)

if(ZLIB_FOUND)
  set_target_properties(tng_io PROPERTIES COMPILE_DEFINITIONS USE_ZLIB=1)
  target_link_libraries(tng_io ${ZLIB_LIBRARIES})
endif()
