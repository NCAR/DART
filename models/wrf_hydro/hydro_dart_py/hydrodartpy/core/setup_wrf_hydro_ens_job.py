import code
import datetime
import os
import pathlib
import subprocess
import sys
import wrfhydropy

def setup_wrf_hydro_ens_job(config, wrf_hydro_ens_sim):

    config_job_exe = config['run_experiment']['job_execution']

    # I think the exe name is hard-coded in wrfhydropy to this. If it's
    # tracked, could bring it in.
    hostname = subprocess.run('hostname', stdout=subprocess.PIPE).stdout.rstrip().decode('utf-8')
    fmt_dict = {
        **{'cmd': './wrf_hydro.exe',
           'hostname': hostname,
           'nproc': config_job_exe['wrf_hydro']['nproc']
        }
    }
    exe_cmd = config_job_exe['wrf_hydro']['exe_cmd'].format(**fmt_dict)

    # This is ugly, but there are no files on disk to query, which is the other way of doing this.
    def get_nlist_date(member):
        file_name_hydro = member.base_hydro_namelist['hydro_nlist']['restart_file']
        date_str_hydro = pathlib.Path(file_name_hydro).name
        date_hydro = datetime.datetime.strptime(date_str_hydro, 'HYDRO_RST.%Y-%m-%d_%H:%M_DOMAIN1')
        if member.base_hrldas_namelist['wrf_hydro_offline']['forc_typ'] not in [9, 10]:
            file_name_hrldas = \
                member.base_hrldas_namelist['noahlsm_offline']['restart_filename_requested']
            date_str_hrldas = pathlib.Path(file_name_hrldas).name
            date_hrldas = datetime.datetime.strptime(date_str_hrldas, 'RESTART.%Y%m%d%H_DOMAIN1')
            if date_hrldas != date_hydro:
                raise ValueError("Namelist restart times do not match.") 
        return date_hydro

    member_dates = [get_nlist_date(mm) for mm in wrf_hydro_ens_sim.members]
    if not all([mm == member_dates[0] for mm in member_dates]):
        raise ValueError("Ensemble members are not at the same times")

    start_time = member_dates[0]
    end_time = \
        member_dates[0] + \
        datetime.timedelta(hours=config['run_experiment']['time']['advance_model_hours'])

    entry_cmd = config['run_experiment']['perturb_forcing']['noise_cmd']
    exit_cmd = None

    job = wrfhydropy.Job(
        exe_cmd=exe_cmd,
        job_id=config_job_exe['scheduler']['job_name'],
        model_start_time=start_time,
        model_end_time=end_time,
        restart=True,
        entry_cmd=entry_cmd,
        exit_cmd=exit_cmd
    )

    wrf_hydro_ens_sim.add(job)

    return wrf_hydro_ens_sim
