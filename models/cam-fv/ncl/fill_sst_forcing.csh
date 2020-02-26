#!/bin/csh 

#PBS  -N fill_sst_forcing
#PBS  -A NCIS0006
#PBS  -q share
# Resources I want:
#    select=#nodes
#    ncpus=#CPUs/node
#    mpiprocs=#MPI_tasks/node
#PBS  -l select=1:ncpus=1:mpiprocs=1
#PBS  -l walltime=03:00:00
# Send email after a(bort) or e(nd)
#PBS  -m ae
#PBS  -M raeder@ucar.edu
# Send standard output and error to this file.
# It's helpful to use the $casename here.
#PBS  -o fill_sst_forcing.eo
#PBS  -j oe 
#--------------------------------------------

# If not in my environment already:
module load gnu
module load ncl

set echo verbose

# set years = ( 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 \
#                                        2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 \
#                                        2010 2011 2012 2013 )
# set years = ( 2018 2019 )
set years = ( 2019 )
set product = avhrr
set dir_out = /glade/scratch/raeder/SST/${product}

cd ${dir_out}
foreach year ($years)
   setenv CYEAR $year
   ncl /gpfs/fs1/work/raeder/Models/CAM_init/SST/Make_yearly_from_daily/fill_sst_forcing.ncl |& tee ncl_${year}.eo
#    ncl fill_sst_forcing.ncl
end

exit
