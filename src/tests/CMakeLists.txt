include_directories(${TRAJECTORY_SOURCE_DIR}/src/lib)
link_directories(${TRAJECTORY_BINARY_DIR}/src/lib)

add_executable(tng_testing tng_io_testing.c)
target_link_libraries(tng_testing tng_io)