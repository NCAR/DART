#!/bin/sh

let test=1

while true; do

 rm -fr testdir
 mkdir testdir
 cd testdir
 ln -s ../file_utils_test .

 case $test in
   1 ) ;;

   2 ) touch dart_log.out  
       chmod 0 dart_log.out
       ;;

   3 ) touch input.nml
       ;;

   4 ) cp ../input.nml .
       chmod 0 input.nml
       ;;

   5 ) echo \&bob > input.nml
       ;;

   6 ) echo \&utilities_nml > input.nml
       ;;

   7 ) echo \&utilities_nml >  input.nml
       echo /               >> input.nml
       ;;

   8 ) echo \&utilities_nml       >  input.nml
       echo /                     >> input.nml
       echo \&file_utils_test_nml >> input.nml
       echo /                     >> input.nml
       ;;

   9 ) echo \&utilities_nml       >  input.nml
       echo \&file_utils_test_nml >> input.nml
       echo /                     >> input.nml
       ;;

  10 ) echo \&utilities_nml       >  input.nml
       echo /                     >> input.nml
       echo \&file_utils_test_nml >> input.nml
       ;;

   * ) exit 
       ;;
 esac
  
cat input.nml
 echo test $test 
 ./file_utils_test

 let test=test+1

 cd ..

done

# only when this is working, clean up
#rmdir -fr testdir

echo done
exit 0


#----------------------------------------------------------------------
# mkdir test; cd test; then:
#  touch dart_log.out; chmod 0 dart_log.out
#  chmod 0 input.nml
#  rm input.nml (no nml file)
#  rm input.nml; touch input.nml (0 length nml file)
#  rm input.nml; echo &bob > input.nml  (no trailing / )
#  rm input.nml; echo &utilities_nml > input.nml  (no trailing / )
#  rm input.nml; echo &file_utils_test_nml > input.nml
#                echo / > input.nml (no &utilities_nml )
# and run this test
#----------------------------------------------------------------------
