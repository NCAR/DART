#!/bin/tcsh

# This is the least inappropriate place for this script.  It uses no matlab,
# but generates files for input to the CAM-SE matlab processor, PlotCubedSpherePatches_eps.m.

# A script to generate the difference between CAM-SE or DART output and the cam_Truth.nc 
# from a perfect model run.
# It should be run in a directory where the output from CAM-SE+DART is living:
#    $case/archive/rest/$date   for CAM-SE initial files
#    $case/archive/dart/hist    for P{oste}rior_Diag.nc files
#

set dates = ( 2005-08-20-00000 )

set truth_dir = /glade/scratch/raeder/SE30_Og16_pmo2/archive/dart/hist
# Set the director, filename root (no date or .nc), and copy (fortran indexing) 
# to extract of the CAM+DART output.
set forec_dir = /glade/scratch/raeder/SE30_Og16_osse5/archive/dart/hist
set forec_file_root = cam_Prior_Diag
set copy = 3

# Make a subdirectory for the intermediate and final output of this script.
if ( $cwd != $forec_dir) cd $forec_dir

foreach date ($dates)
   set f_root = ${forec_file_root}.${date}
   if ( ! -d ${f_root} ) mkdir  ${f_root}
   cd  ${f_root}

   # Exclude the namelist variable from both files because they usually have different sizes.
   ncks -O -F -x -v inputnml -d copy,${copy}  -o ${f_root}.c${copy}_no_nml.nc  ../${f_root}.nc
   echo "Was c${copy}_no_nml created?"
   ls ${f_root}.c${copy}_no_nml.nc
   if (! -f cam_True_State.${date}.no_nml.nc) then
      ncks -O -F -x -v inputnml -o cam_True_State.${date}.no_nml.nc \
                      ${truth_dir}/cam_True_State.${date}.nc
   else
      echo "cam_True_State.${date}.no_nml.nc already exists; using it"
   endif
   
   # Diff the 2 local files.
   ncdiff ${f_root}.c${copy}_no_nml.nc cam_True_State.${date}.no_nml.nc \
          ${forec_file_root}_err.${date}.c${copy}.nc
   echo "ncdiff status = $status"
   echo "Was file created?"
   ls -l ${forec_file_root}_err.${date}.c${copy}.nc

   # Save the True_State file for additional differences.
   if ($status == 0) rm ${f_root}.c${copy}_no_nml.nc
   
   
   echo "Output is in  ${forec_file_root}.${date}, which may not be where you are"
   echo "Run matlab PlotCubedSpherePatches.m, with ${forec_file_root}_err.${date}.c${copy}.nc"
   echo "as the input file, to see (lon,lat) plots of these difference fields."
   
   cd ..
end

exit
