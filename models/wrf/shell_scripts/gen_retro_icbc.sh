#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# DART / WRF gen_retro_icbc – Bash version for PBS (Derecho)
#
# Original behavior:
#   - Loops over times from datea to datefnl
#   - Runs WPS (geogrid/ungrib/metgrid)
#   - Runs real.exe twice (first to make wrfinput+wrfbdy, then to make second-time wrfinput)
#   - Writes:
#       output/${date}/wrfbdy_d01_{gdayf}_{gsecf}_mean
#       output/${date}/wrfinput_d01_{gday}_{gsec}_mean
#       output/${date}/wrfinput_d02_{gday}_{gsec}_mean




#==================================================================
#BSUB -J gen_retro_icbc
#BSUB -o gen_retro_icbc.%J.log
#BSUB -P 25000077
#BSUB -W 0:38
#BSUB -q regular
#BSUB -n 64
#BSUB -x
#BSUB -R "span[ptile=64]"

#PBS -N gen_retro_icbc
#PBS -A 25000077
#PBS -l walltime=00:38:00
#PBS -q regular
#PBS -o gen_retro_icbc.out
#PBS -j oe
#PBS -k eod
#PBS -l select=5:ncpus=60:mpiprocs=60
#PBS -V

set -uo pipefail

echo "gen_retro_icbc.sh is running in $(pwd)"

###############################################################################
# User-configurable dates and param file
###############################################################################

datea=2024051812         # initial cycle time (YYYYMMDDHH)
datefnl=2024052018       # final cycle time   (YYYYMMDDHH)
paramfile="/glade/derecho/scratch/bmraczka/WRFv4.5_nested_bash/scripts/param.sh"

echo "Sourcing parameter file: $paramfile"
source "$paramfile"

# Backups if REMOVE/MOVE/LINK not set in param.sh
: "${REMOVE:=rm -rf}"
: "${MOVE:=mv -f}"
: "${COPY:=cp -p}"
: "${LINK:=ln -sf}"

DART_ADVANCE_TIME="${DART_DIR}/models/wrf/work/advance_time"

###############################################################################
# One-time ICBC directory prep
###############################################################################

cd "$ICBC_DIR"

$REMOVE geo_*.nc namelist.wps namelist.input geogrid_done
mkdir -p geogrid
$LINK "${WPS_SRC_DIR}/geogrid/GEOGRID.TBL" "${ICBC_DIR}/geogrid/GEOGRID.TBL"

mkdir -p "${ICBC_DIR}/metgrid"
$LINK "${WPS_SRC_DIR}/metgrid/METGRID.TBL" "${ICBC_DIR}/metgrid/METGRID.TBL"

###############################################################################
# Helper: run real.exe via PBS job that calls your real.sh
###############################################################################

run_real_via_pbs() {
    local pbs_script="$ICBC_DIR/run_real.pbs"

    cat > "$pbs_script" << EOF
#!/bin/bash
#PBS -N run_real
#PBS -A ${COMPUTER_CHARGE_ACCOUNT}
#PBS -l walltime=00:05:00
#PBS -q ${ADVANCE_QUEUE}
#PBS -l job_priority=${ADVANCE_PRIORITY}
#PBS -o run_real.out
#PBS -j oe
#PBS -k eod
#PBS -l select=1:ncpus=4:mpiprocs=4
#PBS -V

cd "$ICBC_DIR"
"${SHELL_SCRIPTS_DIR}/real.sh" "$paramfile"
EOF

    jobid=$(qsub "$pbs_script")
    echo "Submitted real.exe PBS job: $jobid"
 
    # Wait for the marker file from real.sh
    while [[ ! -e "$ICBC_DIR/real_done" ]]; do
        sleep 15
    done
    rm -f "$ICBC_DIR/real_done"

    # Accumulate log
    if [[ -e "$ICBC_DIR/rsl.out.0000" ]]; then
        cat "$ICBC_DIR/rsl.out.0000" >> "$ICBC_DIR/out.real.exe"
    fi
}

###############################################################################
# Main cycle loop
###############################################################################

while :; do
    echo "    "
    echo "Entering gen_retro_icbc.sh for datea = $datea"

    # Ensure output directory exists for this cycle
    mkdir -p "${OUTPUT_DIR}/${datea}"

    cd "$ICBC_DIR"

    # Link DART input.nml
    $LINK "${RUN_DIR}/input.nml" input.nml

    # Clean old GRIB files
    $REMOVE gfs*pgrb2* *grib2 || true

    # Compute WPS start and end dates
    start_date=$(echo "$datea 0 -w" | "$DART_ADVANCE_TIME")
    end_date=$(echo "$datea 6 -w" | "$DART_ADVANCE_TIME")
    echo "start_date = $start_date"
    echo "end_date   = $end_date"

    # Build namelist.wps via sed script
    cat > script.sed << EOF
/start_date/c\
 start_date = ${start_date},${start_date}
/end_date/c\
 end_date   = ${end_date},${end_date}
EOF

    if [[ "$GRIB_SRC" != "GFS" ]]; then
        echo "ERROR: GRIB_SRC=$GRIB_SRC not supported in this script (expects GFS)."
        exit 2
    fi

    # GRIB file names (GFS 0.25)
    gribfile_a="${GRIB_DATA_DIR}/gfs.0p25.${datea}.f000.grib2"
    gribfile_b="${GRIB_DATA_DIR}/gfs.0p25.${datea}.f006.grib2"

    if [[ ! -r "$gribfile_a" || ! -r "$gribfile_b" ]]; then
        echo "ERROR: GRIB input files not found:"
        echo "  $gribfile_a"
        echo "  $gribfile_b"
        exit 2
    fi

    $LINK "$gribfile_a" GRIBFILE.AAA
    $LINK "$gribfile_b" GRIBFILE.AAB

    sed -f script.sed "${TEMPLATE_DIR}/namelist.wps.template" > namelist.wps
    $LINK "${WPS_SRC_DIR}/ungrib/Variable_Tables/Vtable.${GRIB_SRC}" Vtable

    # Run geogrid once
    if [[ ! -e "${ICBC_DIR}/geogrid_done" ]]; then
        echo "Executing geogrid.exe"
        "${WPS_SRC_DIR}/geogrid.exe" >& output.geogrid.exe
        touch "${ICBC_DIR}/geogrid_done"
    fi

    echo "Executing ungrib.exe"
    $REMOVE output.ungrib.exe."${GRIB_SRC}" || true
    "${WPS_SRC_DIR}/ungrib.exe" >& output.ungrib.exe."${GRIB_SRC}"

    echo "Executing metgrid.exe"
    $REMOVE output.metgrid.exe || true
    "${WPS_SRC_DIR}/metgrid.exe" >& output.metgrid.exe

    # Compute end of assimilation window
    datef=$(echo "$datea $ASSIM_INT_HOURS" | "$DART_ADVANCE_TIME")
    # Gregorian version for wrfbdy naming
    read -r gdayf gsecf _rest < <(echo "$datef 0 -g" | "$DART_ADVANCE_TIME")
    hh=${datea:8:2}

    ###########################################################################
    # Run real.exe twice: first (datea → datef), then (datef → datef)
    ###########################################################################

    for n in 1 2; do
        echo
        echo "RUNNING REAL, STEP $n"
        echo

        if [[ "$n" -eq 1 ]]; then
            date1="$datea"
            date2="$datef"
            fcst_hours="$ASSIM_INT_HOURS"
        else
            date1="$datef"
            date2="$datef"
            fcst_hours=0
        fi

        yyyy1=${date1:0:4}
        mm1=${date1:4:2}
        dd1=${date1:6:2}
        hh1=${date1:8:2}

        yyyy2=${date2:0:4}
        mm2=${date2:4:2}
        dd2=${date2:6:2}
        hh2=${date2:8:2}

        # Build namelist.input from template namelist.input.meso
        cat > script.sed << EOF
/run_hours/c\
    run_hours                  = ${fcst_hours},
/run_minutes/c\
    run_minutes                = 0,
/run_seconds/c\
    run_seconds                = 0,
/start_year/c\
    start_year                 = ${yyyy1}, ${yyyy1},
/start_month/c\
    start_month                = ${mm1}, ${mm1},
/start_day/c\
    start_day                  = ${dd1}, ${dd1},
/start_hour/c\
    start_hour                 = ${hh1}, ${hh1},
/start_minute/c\
    start_minute               = 00, 00,
/start_second/c\
    start_second               = 00, 00,
/end_year/c\
    end_year                   = ${yyyy2}, ${yyyy2},
/end_month/c\
    end_month                  = ${mm2}, ${mm2},
/end_day/c\
    end_day                    = ${dd2}, ${dd2},
/end_hour/c\
    end_hour                   = ${hh2}, ${hh2},
/end_minute/c\
    end_minute                 = 00, 00,
/end_second/c\
    end_second                 = 00, 00,
EOF

        sed -f script.sed "${TEMPLATE_DIR}/namelist.input.meso" > namelist.input

        # Clean flags & logs
        $REMOVE real_done rsl.* script.sed || true
        : > out.real.exe

        # Submit PBS job that runs real.exe via real.sh
        run_real_via_pbs

        # On return, move wrfinput/wrfbdy to appropriate locations
        read -r gday gsec _rest < <(echo "$date1 0 -g" | "$DART_ADVANCE_TIME")

        # Move wrfinput_d01/d02
        if [[ -e wrfinput_d01 ]]; then
            $MOVE wrfinput_d01 "${OUTPUT_DIR}/${datea}/wrfinput_d01_${gday}_${gsec}_mean"
        fi
        if [[ -e wrfinput_d02 ]]; then
            $MOVE wrfinput_d02 "${OUTPUT_DIR}/${datea}/wrfinput_d02_${gday}_${gsec}_mean"
        fi

        # For first real.exe call, also move wrfbdy_d01
        if [[ "$n" -eq 1 && -e wrfbdy_d01 ]]; then
            $MOVE wrfbdy_d01 "${OUTPUT_DIR}/${datea}/wrfbdy_d01_${gdayf}_${gsecf}_mean"
        fi

    done  # end n=1,2

    ###########################################################################
    # Move to next time, or exit if final time reached
    ###########################################################################
    if [[ "$datea" == "$datefnl" ]]; then
        echo "Reached final date $datefnl – script exiting normally."
        exit 0
    fi

    # Advance datea by ASSIM_INT_HOURS
    datea=$(echo "$datea $ASSIM_INT_HOURS" | "$DART_ADVANCE_TIME" | awk '{print $1}')
    echo "     "
    echo "Starting next time: $datea"

done

exit 0

