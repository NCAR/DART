#!/bin/bash
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
# AUTHOR:   T. R. Whitcomb
#           Naval Research Laboratory
#
# Allows the use of files in DART that are not included with the standard DART
# distribution.  Use of these files should be minimized, as they are fairly
# inflexible and require large changes whenever DART is upgraded.
#
# This script looks for directories in the "externals" subdirectory that
# mirror the DART directory structure of where they actually belong.  It moves
# the original file out of the way if necessary, then installs a symbolic link
# to the updated file.
#
# The command line argument supplies the location of the DART root directory,
# which is used to find all the files.
######

######
# BEGIN DEFINITIONS
######

EXTERNALS_PATH=models/coamps/externals

# process_external takes a file name as an input, gets the target information,
# then creates a link to the new file in the proper directory
process_external()
{
    modified_file="${1}"
    if [[ -z "${modified_file}" ]]
    then
        echo "No filename supplied!"
        exit 1
    fi

    modified_dir="${2}"
    if [[ -z "${modified_dir}" ]]
    then    
        echo "No directory supplied!"
        exit 1
    fi

    # No need to process README files
    if [[ "${modified_file}" == "README" ]]
    then 
        echo "Skipping README..."
        return 0
    fi

    # Since the externals directory structure mirrors the main
    # DART directory structure, to get the target file path we
    # just need to delete the bit that points to the externals
    # e.g. DART/models/coamps/externals/ -> DART/   
    target_dir="${modified_dir/${EXTERNALS_PATH}/}"

    # This is necessary since we may be adding directories that
    # do not exist in the basic DART distribution (this is mainly
    # for location modules)
    mkdir -p "${target_dir}"

    target_file="${modified_file}"
    target_file_fullpath="${target_dir}/${target_file}"
    modified_file_fullpath="${modified_dir}/${modified_file}"

    # Just save a hidden copy instead of overwriting or changing
    # the name to anything fancy, but just delete it if it's already
    # a symbolic link
    if [[ -e "${target_file_fullpath}" ]] 
    then
        if [[ -h "${target_file_fullpath}" ]]
        then
            echo "Deleting symbolic link ${target_file}..."
            rm "${target_file_fullpath}"
        else
            echo "Archiving ${target_file}..."
            archive_file="${target_dir}/.${target_file}"
            mv "${target_file_fullpath}" "${archive_file}" 
        fi
    fi  

    # Now actually do the linking
    ln -s "${modified_file_fullpath}" "${target_file_fullpath}"
    echo " Linked in ${modified_file} to ${target_file_fullpath}"
}

# traverse_tree takes a root directory as an input, then recursively descends
# through the directory tree 
traverse_tree()
{
    target_dir="${1:-.}"

    pushd "$target_dir" > /dev/null
    current_dir=`pwd`
    for to_process in *
    do
        if [[ -f "${to_process}" ]]
        then
            # Process files
            process_external "${to_process}" "${current_dir}"
        elif [[ -d "${to_process}" ]]
        then
            # Descend into directories
            traverse_tree "${to_process}" "${current_dir}"
        fi
    done
    popd > /dev/null
}

######
# END DEFINITIONS
######


# If the current directory isn't supplied, assume that we're running from the
# shell_scripts directory
supplied_dart_dir="$1"
calculated_dart_dir=`cd ../../../;pwd`
dart_dir="${supplied_dart_dir:-${calculated_dart_dir}}"
echo "Using DART root directory $dart_dir"

# Start processing the externals - as long as the directory structure in
# "externals" mirrors that of the main DART directory, this should work
traverse_tree "${dart_dir}/${EXTERNALS_PATH}"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

