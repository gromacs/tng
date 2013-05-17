set(COMP_INC_DIR ${CMAKE_SOURCE_DIR}/include/compression)

include_directories(${COMP_INC_DIR})

set(HEADER_FILES ${COMP_INC_DIR}/bwlzh.h ${COMP_INC_DIR}/bwt.h ${COMP_INC_DIR}/coder.h ${COMP_INC_DIR}/dict.h ${COMP_INC_DIR}/fixpoint.h ${COMP_INC_DIR}/huffman.h ${COMP_INC_DIR}/lz77.h ${COMP_INC_DIR}/merge_sort.h ${COMP_INC_DIR}/mtf.h ${COMP_INC_DIR}/my64bit.h ${COMP_INC_DIR}/rle.h ${COMP_INC_DIR}/tng_compress.h ${COMP_INC_DIR}/vals16.h ${COMP_INC_DIR}/warnmalloc.h ${COMP_INC_DIR}/widemuldiv.h)

add_library(tng_compress SHARED bwlzh.c bwt.c coder.c dict.c fixpoint.c huffman.c huffmem.c lz77.c merge_sort.c mtf.c rle.c tng_compress.c vals16.c warnmalloc.c widemuldiv.c xtc2.c xtc3.c)
if(UNIX)
target_link_libraries(tng_compress m)
endif()
