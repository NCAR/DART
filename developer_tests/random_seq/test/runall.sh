#!/bin/sh

for i in test_*
do
  echo running $i
  ./$i
done

exit 0
