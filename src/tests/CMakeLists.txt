include_directories(${TRAJECTORY_SOURCE_DIR}/src/lib)
link_directories(${TRAJECTORY_BINARY_DIR}/src/lib)

add_executable(trg_testing trg_io_testing.c)
target_link_libraries(trg_testing trg_io)