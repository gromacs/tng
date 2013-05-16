link_directories(${TRAJECTORY_BINARY_DIR}/src/lib)

set(TNG_COMPRESS_FILES_DIR ${CMAKE_BINARY_DIR}/test_tng_compress_files/ CACHE STRING "Directory where to write tng_compress test files.")

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/test_tng_compress_files)
configure_file(${TRAJECTORY_SOURCE_DIR}/include/compression/tng_compress_testing.h.in ${CMAKE_BINARY_DIR}/generated/tng_compress_testing.h)

include_directories(${CMAKE_SOURCE_DIR}/include/compression)
include_directories(${CMAKE_BINARY_DIR}/generated/)

set(number 0)
set(numtests 62)

while( number LESS ${numtests})

math( EXPR number "${number} + 1" )

add_executable(test_tng_compress_gen${number} testsuite.c)
target_link_libraries(test_tng_compress_gen${number} tng_compress m)
list(APPEND gen${number}_build_definitions GEN)
list(APPEND gen${number}_build_definitions TESTPARAM="test${number}.h")
set_target_properties(test_tng_compress_gen${number} PROPERTIES COMPILE_DEFINITIONS "${gen${number}_build_definitions}")
add_dependencies(test_tng_compress_gen${number} test${number}.h)

add_executable(test_tng_compress_read${number} testsuite.c)
target_link_libraries(test_tng_compress_read${number} tng_compress m)
list(APPEND read${number}_build_definitions TESTPARAM="test${number}.h")
set_target_properties(test_tng_compress_read${number} PROPERTIES COMPILE_DEFINITIONS "${read${number}_build_definitions}")
add_dependencies(test_tng_compress_read${number} test${number}.h)

endwhile()

file(COPY test_tng_compress_write.sh DESTINATION ${CMAKE_BINARY_DIR}/bin)
file(COPY test_tng_compress_read.sh DESTINATION ${CMAKE_BINARY_DIR}/bin)
file(COPY test_tng_compress_write.bat DESTINATION ${CMAKE_BINARY_DIR}/bin)
file(COPY test_tng_compress_read.bat DESTINATION ${CMAKE_BINARY_DIR}/bin)
