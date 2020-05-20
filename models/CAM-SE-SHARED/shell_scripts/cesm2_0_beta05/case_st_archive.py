# This code may (or may not) be part of the CESM distribution,
# So it is not protected by the DART copyright agreement. 
# This file is a variation of 
# cesm2_0_beta05/cime/utils/python/CIME/case_st_archive.py

# DART $Id$

"""
short term archiving
"""

import shutil, glob, re, os

from CIME.XML.standard_module_setup import *
from CIME.case_submit               import submit
from CIME.XML.env_archive           import EnvArchive
from CIME.utils                     import append_status
from os.path                        import isdir, join

logger = logging.getLogger(__name__)

###############################################################################
def _get_datenames(case):
###############################################################################

    logger.debug('In get_datename...')
    rundir = case.get_value('RUNDIR')
    expect(isdir(rundir), 'Cannot open directory %s ' % rundir)
    casename = case.get_value("CASE")
    files = sorted(glob.glob(os.path.join(rundir, casename + '.cpl.r*.nc')))
    if not files:
        expect(False, 'Cannot find a %s.cpl.r.*.nc file in directory %s ' % (casename, rundir))
    datenames = []
    for filename in files:
        names = filename.split('.')
        datename = names[-2]
        datenames.append(datename)
        logger.debug('cpl dateName: %s ' % datename)
    return datenames


###############################################################################
def _get_ninst_info(case, compclass):
###############################################################################

    if compclass != 'cpl':
        ninst = case.get_value('NINST_' + compclass.upper())
    else:
        ninst = 1
    ninst_strings = []
    if ninst is None:
        ninst = 1
    for i in range(1,ninst+1):
        if ninst > 1:
           ninst_strings.append('_' + '%04d' % i)
# KDR This doesn't do what's intended.  The result is ninst_strings is not empty.
#     It has "''" in it.  Leave it for now, since _archive_restarts does what I need it too.
##        else:
##            ninst_strings.append('')

    logger.debug("ninst and ninst_strings are: %s and %s for %s" %(ninst, ninst_strings, compclass))
    return ninst, ninst_strings


###############################################################################
def _archive_rpointer_files(case, archive, archive_entry, archive_restdir,
                            datename, datename_is_last):
###############################################################################

    # archive the rpointer files associated with datename
    casename = case.get_value("CASE")
    compclass = archive.get_entry_info(archive_entry)[1]
    ninst_strings = _get_ninst_info(case, compclass)[1]

    if datename_is_last:
        # Copy of all rpointer files for latest restart date
        rundir = case.get_value("RUNDIR")
        rpointers = glob.glob(os.path.join(rundir, 'rpointer.*'))
        for rpointer in rpointers:
            shutil.copy(rpointer, os.path.join(archive_restdir, os.path.basename(rpointer)))
    else:
        # Generate rpointer file(s) for interim restarts for the one datename and each
        # possible value of ninst_strings
        if case.get_value('DOUT_S_SAVE_INTERIM_RESTART_FILES'):

            # parse env_archive.xml to determine the rpointer files
            # and contents for the given archive_entry tag
            rpointer_items = archive.get_rpointer_contents(archive_entry)

            # loop through the possible rpointer files and contents
            for rpointer_file, rpointer_content in rpointer_items:

                # put in a temporary setting for ninst_strings if they are empty
                # in order to have just one loop over ninst_strings below
                if rpointer_content is not 'unset':
                    if not ninst_strings:
                        ninst_strings = ["empty"]
                    for ninst_string in ninst_strings:
                        if ninst_string == 'empty':
                            ninst_string = ""
                        for key, value in [('$CASE', casename),
                                           ('$DATENAME', datename),
                                           ('$NINST_STRING', ninst_string)]:
                            rpointer_file = rpointer_file.replace(key, value)
                            rpointer_content = rpointer_content.replace(key, value)
    
                        # write out the respect files with the correct contents
                        rpointer_file = os.path.join(archive_restdir, rpointer_file)
                        logger.info("writing rpointer_file %s" % rpointer_file)
                        f = open(rpointer_file, 'w')
                        for output in rpointer_content.split(','):
                            f.write("%s \n" %output)
                        f.close()


###############################################################################
def _archive_log_files(case):
###############################################################################

    dout_s_root = case.get_value("DOUT_S_ROOT")
    rundir = case.get_value("RUNDIR")
    archive_logdir = os.path.join(dout_s_root, 'logs')
    if not os.path.exists(archive_logdir):
        os.makedirs(archive_logdir)
        logger.debug("created directory %s " %archive_logdir)

# KDR want dart_logs too
#    logfiles = glob.glob(os.path.join(rundir, '*.log.*'))
    logfiles = glob.glob(os.path.join(rundir, '*[_.]log.*'))
    logger.debug("short term archiving logfiles are")
    logger.debug("  %s" %(logfiles))

    for logfile in logfiles:
        srcfile = join(rundir, os.path.basename(logfile))
        destfile = join(archive_logdir, os.path.basename(logfile))
        shutil.move(srcfile, destfile)
        logger.info("moving \b%s to \b%s" %(srcfile, destfile))


###############################################################################
def _archive_history_files(case, archive, archive_entry,
                           compclass, compname, histfiles_savein_rundir):
###############################################################################
    """
    perform short term archiving on history files in rundir
    """

    # determine history archive directory (create if it does not exist)
    dout_s_root = case.get_value("DOUT_S_ROOT")
    casename = case.get_value("CASE")
    archive_histdir = os.path.join(dout_s_root, compclass, 'hist')
    if not os.path.exists(archive_histdir):
        os.makedirs(archive_histdir)
        logger.debug("created directory %s" %archive_histdir)

    # determine ninst and ninst_string
    ninst, ninst_string = _get_ninst_info(case, compclass)

    # archive history files - the only history files that kept in the
    # run directory are those that are needed for restarts
    rundir = case.get_value("RUNDIR")
    for suffix in archive.get_hist_file_extensions(archive_entry):
        for i in range(ninst):
# KDR I've chosen (so far) to  name the (non-NetCDF) obs_seq files 
#     with no $case.$comp at their beginnings.
#     Modify this section to make it see those files.
#     I also chose to use compclass 'cam' for DART output from CAM assims
#       to distinguish it from other components in multi-component assims.
#             if compname == 'dart':
            if compclass == 'esp' and compname == 'cam':
#               newsuffix = casename + suffix
                newsuffix = casename + '.' + compname + ".*" + suffix
                if re.compile('_obs_seq|dart_log').search(suffix):
                   newsuffix = suffix
# KDR end of new
            else:
                if ninst_string:
                    newsuffix = casename + '.' + compname + ".*" + ninst_string[i] + suffix
                else:
                    newsuffix = casename + '.' + compname + ".*" + suffix
            logger.debug("short term archiving suffix is %s " %newsuffix)
            pfile = re.compile(newsuffix)
            histfiles = [f for f in os.listdir(rundir) if pfile.search(f)]
            if histfiles:
                logger.debug("hist files are %s " %histfiles)
                for histfile in histfiles:
                    srcfile = join(rundir, histfile)
                    expect(os.path.isfile(srcfile),
                           "history file %s does not exist " %srcfile)
                    destfile = join(archive_histdir, histfile)
                    if histfile in histfiles_savein_rundir:
                        logger.info("histfiles_savein_rudir copying \n%s to \n%s\n " %(srcfile, destfile))
                        shutil.copy(srcfile, destfile)
                    else:
                        logger.info("histfiles moving \n%s to \n%s " %(srcfile, destfile))
                        shutil.move(srcfile, destfile)


###############################################################################
def get_histfiles_for_restarts(case, archive, archive_entry, restfile):
###############################################################################

    # determine history files that are needed for restarts
    histfiles = []
    rest_hist_varname = archive.get_entry_value('rest_history_varname', archive_entry)
    if rest_hist_varname != 'unset':
        rundir = case.get_value("RUNDIR")
        cmd = "ncdump -v %s %s " %(rest_hist_varname, os.path.join(rundir, restfile))
        rc, out, error = run_cmd(cmd)
        if rc != 0:
            logger.debug(" WARNING: %s failed rc=%d\nout=%s\nerr=%s" %(cmd, rc, out, error))

        searchname = "%s =" %rest_hist_varname
        if searchname in out:
            offset = out.index(searchname)
            items = out[offset:].split(",")
            for item in items:
                # the following match has an option of having a './' at the beginning of
                # the history filename
                matchobj = re.search(r"\"(\.*\/*\w.*)\s?\"", item)
                if matchobj:
                    histfile = matchobj.group(1).strip()
                    histfile = os.path.basename(histfile)
                    # append histfile to the list ONLY if it exists in rundir before the archviving
                    if os.path.isfile(os.path.join(rundir,histfile)): 
                        histfiles.append(histfile)
    return histfiles


###############################################################################
def _archive_restarts(case, archive, archive_entry,
                      compclass, compname, datename, datename_is_last):
###############################################################################

    # determine directory for archiving restarts based on datename
    dout_s_root = case.get_value("DOUT_S_ROOT")
    rundir = case.get_value("RUNDIR")
    casename = case.get_value("CASE")
    archive_restdir = join(dout_s_root, 'rest', datename)
    if not os.path.exists(archive_restdir):
        os.makedirs(archive_restdir)

    # archive the rpointer file(s) for this datename and all possible ninst_strings
    _archive_rpointer_files(case, archive, archive_entry, archive_restdir,
                            datename, datename_is_last)

    # determine ninst and ninst_string
    ninst, ninst_strings = _get_ninst_info(case, compclass)

    # move all but latest restart files into the archive restart directory
    # copy latest restart files to archive restart directory
    histfiles_savein_rundir = []

    # get file_extension suffixes
    for suffix in archive.get_rest_file_extensions(archive_entry):
        for i in range(ninst):
            restfiles = ""
            pattern = r"%s\.%s\d*.*" % (casename, compname)
# KDR Looking for DART
#             logger.info("Creating rest_files; pattern, suffix, compname: %s %s %s" %(pattern, suffix, compname))

            if compclass == 'esp' and 'cam' in pattern:
                pattern = suffix + datename
# KDR Print what's being archived for DART
#                 logger.info("Archiving DART suffix %s" %(suffix))
                pfile = re.compile(pattern)
                restfiles = [f for f in os.listdir(rundir) if pfile.search(f)]
            else:
                pfile = re.compile(pattern)
                files = [f for f in os.listdir(rundir) if pfile.search(f)]
                if ninst_strings:
                    pattern = ninst_strings[i] + suffix + datename
                else:
                    pattern = suffix + datename

                logger.info("restfiles for patt,inst,suff,date %s %s %s" %(pattern,  suffix, datename))
                logger.info("    %s\n" %(restfiles))
                pfile = re.compile(pattern)
                restfiles = [f for f in files if pfile.search(f)]

            for restfile in restfiles:
                restfile = os.path.basename(restfile)

                # obtain array of history files for restarts
                # need to do this before archiving restart files
                histfiles_for_restart = get_histfiles_for_restarts(case, archive,
                                                                   archive_entry, restfile)

# KDR Looking for extra copys
#                 logger.info("histfiles_for_restart for %s %s %s" %(pattern, suffix, datename))
#                 logger.info("    %s\n" %(histfiles_for_restart))

                if datename_is_last and histfiles_for_restart:
                    for histfile in histfiles_for_restart:
                        if histfile not in histfiles_savein_rundir:
                            histfiles_savein_rundir.append(histfile)

                # archive restart files and all history files that are needed for restart
                # Note that the latest file should be copied and not moved
                if datename_is_last:
                    srcfile = os.path.join(rundir, restfile)
                    destfile = os.path.join(archive_restdir, restfile)
                    shutil.copy(srcfile, destfile)
                    logger.info("datename_is_last copying \n%s to \n%s\n" %(srcfile, destfile))
                    for histfile in histfiles_for_restart:
                        srcfile = os.path.join(rundir, histfile)
                        destfile = os.path.join(archive_restdir, histfile)
                        expect(os.path.isfile(srcfile),
                               "restart file %s does not exist " %srcfile)
                        shutil.copy(srcfile, destfile)
                        logger.info("datename_is_last + histfiles_for_restart copying \n%s to \n%s\n" %(srcfile, destfile))
                else:
                    # Only archive intermediate restarts if requested - otherwise remove them
                    if case.get_value('DOUT_S_SAVE_INTERIM_RESTART_FILES'):
                        srcfile = os.path.join(rundir, restfile)
                        destfile = os.path.join(archive_restdir, restfile)
# KDR redundant?                        logger.info("moving \n%s to \n%s" %(srcfile, destfile))
                        expect(os.path.isfile(srcfile),
                               "restart file %s does not exist " %srcfile)
                        shutil.move(srcfile, destfile)
                        logger.info("save_interim_restarts moving \n%s to \n%s" %(srcfile, destfile))

                        # need to copy the history files needed for interim restarts - since
                        # have not archived all of the history files yet
                        for histfile in histfiles_for_restart:
                            srcfile = os.path.join(rundir, histfile)
                            destfile = os.path.join(archive_restdir, histfile)
                            expect(os.path.isfile(srcfile),
                                   "hist file %s does not exist " %srcfile)
                            shutil.copy(srcfile, destfile)
                            logger.info("histfiles_for_restart copying \n%s to \n%s\n" %(srcfile, destfile))
                    else:
                        srcfile = os.path.join(rundir, restfile)
                        logger.info("removing interim restart file %s" %srcfile)
                        if (os.path.isfile(srcfile)):
                            try:
                                os.remove(srcfile)
                            except OSError:
                                logger.warn("unable to remove interim restart file %s" %srcfile)
                        else:
                            logger.warn("interim restart file %s does not exist" %srcfile)

    return histfiles_savein_rundir

###############################################################################
def _archive_process(case, archive):
###############################################################################
    """
    Parse config_archive.xml and perform short term archiving
    """

    logger.debug('In archive_process...')
    compset_comps = case.get_compset_components()
    compset_comps.append('cpl')
    compset_comps.append('dart')

    # archive log files
    _archive_log_files(case)

    for archive_entry in archive.get_entries():
        # determine compname and compclass
        compname, compclass = archive.get_entry_info(archive_entry)

        # check for validity of compname
        if compname not in compset_comps:
# KDR Say if something is being skipped
            logger.info('Skipping compname %s' %(compname))
            continue

        # archive restarts and all necessary associated fields (e.g. rpointer files)
        logger.info('-------------------------------------------')
        logger.info('doing short term archiving for %s (%s)' % (compname, compclass))
        logger.info('-------------------------------------------')
        datenames = _get_datenames(case)
        for datename in datenames:
            datename_is_last = False
            if datename == datenames[-1]:
                datename_is_last = True

            # archive restarts
            histfiles_savein_rundir = _archive_restarts(case, archive, archive_entry,
                                                        compclass, compname,
                                                        datename, datename_is_last)

            # if the last datename for restart files, then archive history files
            # for this compname
            if datename_is_last:
                logger.info("histfiles_savein_rundir %s " %histfiles_savein_rundir)
                _archive_history_files(case, archive, archive_entry,
                                       compclass, compname, histfiles_savein_rundir)

###############################################################################
def restore_from_archive(case):
###############################################################################
    """
    # NOT CLEAR IF THIS IS CALLED YET
    Take most recent archived restart files and load them into current case.
    """
    dout_sr = case.get_value("DOUT_S_ROOT")
    rundir = case.get_value("RUNDIR")
    most_recent_rest = run_cmd_no_fail("ls -1dt %s/rest/* | head -1" % dout_sr)

    for item in glob.glob("%s/*" % most_recent_rest):
        base = os.path.basename(item)
        dst = os.path.join(rundir, base)
        if os.path.exists(dst):
            os.remove(dst)

        shutil.copy(item, rundir)


###############################################################################
def case_st_archive(case):
###############################################################################
    """
    Create archive object and perform short term archiving
    """
    caseroot = case.get_value("CASEROOT")

    dout_s_root = case.get_value('DOUT_S_ROOT')
    if dout_s_root is None or dout_s_root == 'UNSET':
        expect(False,
               'XML variable DOUT_S_ROOT is required for short-term achiver')
    if not isdir(dout_s_root):
        os.makedirs(dout_s_root)

    dout_s_save_interim = case.get_value('DOUT_S_SAVE_INTERIM_RESTART_FILES')
    if dout_s_save_interim == 'FALSE' or dout_s_save_interim == 'UNSET':
        rest_n = case.get_value('REST_N')
        stop_n = case.get_value('STOP_N')
        if rest_n < stop_n:
            logger.warn('Restart files from end of run will be saved'
                        'interim restart files will be deleted')

    logger.info("st_archive starting")

    # do short-term archiving
    append_status("st_archiving starting", caseroot=caseroot, sfile="CaseStatus")

    archive = EnvArchive(infile=os.path.join(caseroot, 'env_archive.xml'))

    _archive_process(case, archive)

    append_status("st_archiving completed", caseroot=caseroot, sfile="CaseStatus")
    logger.info("st_archive completed")

    # resubmit case if appropriate
    resubmit = case.get_value("RESUBMIT")
    if resubmit > 0:
        append_status("resubmitting from st_archive",
                      caseroot=caseroot, sfile="CaseStatus")
        logger.info("resubmitting from st_archive, resubmit=%d"%resubmit)
        if case.get_value("MACH") == "mira":
            expect(os.path.isfile(".original_host"), "ERROR alcf host file not found")
            with open(".original_host", "r") as fd:
                sshhost = fd.read()
            run_cmd("ssh cooleylogin1 ssh %s '%s/case.submit %s --resubmit' "\
                        %(sshhost, caseroot, caseroot), verbose=True)
        else:
            submit(case, resubmit=True)

    return True

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$ -->
