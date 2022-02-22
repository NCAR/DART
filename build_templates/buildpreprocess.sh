#-------------------------
# Build and run preprocess
# Arguements: 
#  none
# Globals:
#  DART - root of DART
#-------------------------
function buildpreprocess() {

 local pp_dir=$DART/assimilation_code/programs/preprocess

 # run preprocess if it is in the current directory
 if [ -f preprocess ]; then
   echo "already there"
   ./preprocess
   return
 fi

# link to preprocess if it is already built
if [ -f $pp_dir/preprocess ]; then
   echo "not there but built"
   ln -s $pp_dir/preprocess .
   ./preprocess 
   return
fi

 # build preproces
 cd $pp_dir
 $DART/build_templates/mkmf -x -p $pp_dir/preprocess \
      -a $DART $pp_dir/path_names_preprocess
 cd -
 ln -s $pp_dir/preprocess .
 ./preprocess
}

