from shutil import copyfile, copyfileobj
import getpass
import os
import sys
import glob
from pathlib import Path

# -----------------------------------------------------#
#                 Define needed function               #
# -----------------------------------------------------#


def unzipNsave(old, new):
    '''
    Saves *.gz file in old as unzipped file in new
    '''
    import gzip

    with gzip.open(old, 'rb') as f_in:
        with open(new, 'wb') as f_out:
            copyfileobj(f_in, f_out)


# -----------------------------------------------------#
#                Main body of the script               #
# -----------------------------------------------------#

#                    Get username                      #
# -----------------------------------------------------#

USER = getpass.getuser()

#      Create stage directory if it doesn't exist      #
# -----------------------------------------------------#

case_name = 'g210.G_JRA.v14.gx1v7.01'
stagedir = '/glade/scratch/' + USER + '/' + case_name + '/rest/2010-01-01-00000/'
Path("stagedir").mkdir(parents=True, exist_ok=True)

#    First load from the directory of MONTHLY saves    #
# -----------------------------------------------------#

load_path = '/glade/campaign/cgd/oce/people/whokim/csm/' + case_name + '/rest/'

years = list(range(41, 202, 10))+list(range(209, 251))
# For debugging:
# years = list(range(31,32))
inds = list(range(1, len(years)+1))

print('Copying from ' + load_path + ' to ' + stagedir)

for ii, year in enumerate(years):

    # Zero pad the year and index string so that there are 4 digits
    zero_filled_year = str(year).zfill(4)
    zero_filled_ii = str(inds[ii]).zfill(4)

    print('Copying year ' + str(year) + ' to ensemble number '
          + zero_filled_ii)

    yearpath = load_path + '/' + zero_filled_year + '-01-01-00000/'

    # Define file paths
    popold = yearpath + case_name + '.pop.r.' + zero_filled_year \
        + '-01-01-00000.nc'
    ovfold = yearpath + case_name + '.pop.ro.' + zero_filled_year \
        + '-01-01-00000'
    ciceold = yearpath + case_name + '.cice.r.' + zero_filled_year \
        + '-01-01-00000.nc'

    popnew = stagedir + case_name + '.pop_' + zero_filled_ii \
        + '.r.2010-01-01-00000.nc'
    ovfnew = stagedir + case_name + '.pop_' + zero_filled_ii \
        + '.ro.2010-01-01-00000'
    cicenew = stagedir + case_name + '.cice_' + zero_filled_ii \
        + '.r.2010-01-01-00000.nc'

    # Copy over
    copyfile(popold, popnew)
    copyfile(ovfold, ovfnew)
    copyfile(ciceold, cicenew)


#         Then load the rest from DAILY saves
#      Some years (251-262) must be decompressed
# -----------------------------------------------------#

case_name_dr = 'g210.G_JRA.v14.gx1v7.01.dr'
load_path = '/glade/campaign/cgd/oce/people/whokim/csm/' + case_name_dr + '/rest'

ly = len(years)
years2 = list(range(251, 272))

# Start from where we left off with the daily saves with inds
inds = list(range(ly+1, ly+len(years2)+1))

print('Copying from ' + load_path + ' to ' + stagedir)

for ii, year in enumerate(years2):

    # Zero pad the year and index string so that there are 4 digits
    zero_filled_year = str(year).zfill(4)
    zero_filled_ii = str(inds[ii]).zfill(4)

    yearpath = load_path + '/' + zero_filled_year + '-01-01-00000/'

    # Define file paths
    popold = glob.glob(yearpath + case_name_dr + '.pop.r.' + zero_filled_year
                       + '-01-01-00000.nc*')[0]
    ovfold = glob.glob(yearpath + case_name_dr + '.pop.ro.' + zero_filled_year
                       + '-01-01-00000*')[0]
    ciceold = glob.glob(yearpath + case_name_dr + '.cice.r.' + zero_filled_year
                        + '-01-01-00000.nc*')[0]

    popnew = stagedir + case_name + '.pop_' + zero_filled_ii \
        + '.r.2010-01-01-00000.nc'
    ovfnew = stagedir + case_name + '.pop_' + zero_filled_ii \
        + '.ro.2010-01-01-00000'
    cicenew = stagedir + case_name + '.cice_' + zero_filled_ii \
        + '.r.2010-01-01-00000.nc'

    # Copy over, decompressing as necessary

    if popold[-3:] == '.gz':
        print('Decompressing and copying year ' + str(year)
              + ' to ensemble number ' + zero_filled_ii)
    else:
        print('Copying year ' + str(year) + ' to ensemble number '
              + zero_filled_ii)

    if popold[-3:] == '.gz':
        unzipNsave(popold, popnew)
    else:
        copyfile(popold, popnew)

    # OVF
    if ovfold[-3:] == '.gz':
        unzipNsave(ovfold, ovfnew)
    else:
        copyfile(ovfold, ovfnew)

    # CICE
    if ciceold[-3:] == '.gz':
        unzipNsave(ciceold, cicenew)
    else:
        copyfile(ciceold, cicenew)

print('Completed')
