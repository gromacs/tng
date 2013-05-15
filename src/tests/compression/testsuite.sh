#!/bin/sh
do_write_test="Yes"
if [ -n "$1" ]; then
   do_write_test=""
fi
STARTTEST=1
ENDTEST=57
#CFLAGS="-Wall -O2 -I../../../include/compress"
CFLAGS="-O2 -I../../../include/compression"
LIBS="-L../../../build/lib -ltng_compress -lm"
LD_LIBRARY_PATH="../../../build/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH
#CFLAGS="-O0 -Wall -g"
#LIBS="-lm -lefence"
CC="gcc"
# 32 bit
#CC="gcc -m32"
for testnum in $(seq $STARTTEST $ENDTEST); do
    testname=$(grep "TESTNAME" test$testnum.h|sed 's/#define TESTNAME//')
    sed "s/TESTPARAM/\"test$testnum.h\"/" <testsuite.c >test$testnum.c
    if [ -n "$do_write_test" ]; then
        echo Write test $testnum: $testname
        $CC -DGEN $CFLAGS -o gen$testnum test$testnum.c -ltng_compress $LIBS
        ./gen$testnum
        rm -f gen$testnum
    fi
    echo Read test $testnum: $testname
    $CC $CFLAGS -o read$testnum test$testnum.c -ltng_compress $LIBS
    ./read$testnum
    rm -f read$testnum
    rm -f test$testnum.c
done

    
    
