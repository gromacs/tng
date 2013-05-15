#!/bin/sh
STARTTEST=1
ENDTEST=57
for testnum in $(seq $STARTTEST $ENDTEST); do
    if [ -r test$testnum.tng ]; then
      grep -v "EXPECTED_FILESIZE" test$testnum.h >tmp$$.h
      echo "#define EXPECTED_FILESIZE" $(ls -l test$testnum.tng |awk '{print $5}'). >>tmp$$.h
      mv tmp$$.h test$testnum.h
    fi
done
