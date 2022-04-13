import copy
from datetime import datetime, timedelta
import f90nml
import math
import os
import pathlib
import pickle
import shlex
import shutil
import subprocess
import uuid
import warnings

from wrfhydropy.core.model import get_git_revision_hash
from wrfhydropy.core.ensemble_tools import get_ens_dotfile_end_datetime

from wrfhydropy import Job, Scheduler


class DartWork(object):
    def __init__(
        self,
        path_rel,
        path_dart,
        build_dir,
        mpi
    ):

        self.path_rel = path_rel

        self.input_nml_file = path_dart / self.path_rel / 'input.nml'
        self.input_nml = f90nml.read(self.input_nml_file)

        if build_dir != path_dart:
            self.input_nml_file = build_dir / self.path_rel / 'input.nml'
            (build_dir / self.path_rel).mkdir(exist_ok=True, parents=True)
            self.input_nml.write(self.input_nml_file)

        self.compile(path_dart, path_rel, mpi)

        # list the mkmfs and get the binaries.
        # HK note: the mkmf_ files no longer exist. 
        mkmfs = list((pathlib.PosixPath(path_dart) / path_rel).glob('mkmf_*'))
        dart_exes = [pathlib.PosixPath(str(mm).replace('/mkmf_','/')) for mm in mkmfs]
        build_exes = [(build_dir / path_rel / dd.name) for dd in dart_exes]
        _ = [shutil.copy(str(dd), str(ss)) for dd, ss in zip(dart_exes, build_exes)] 
        self.exes = {bb.name: bb for bb in build_exes}


    def compile(
        self,
        path_dart,
        path_rel,
        mpi
    ):

        # compile each workdir
        build_cmd = './quickbuild.sh'
        if not mpi:
            build_cmd += 'nompi'

        print(path_rel + ': Running "' + build_cmd + '"')
        self.compile_subproc = subprocess.run(
            shlex.split(build_cmd),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=path_dart / path_rel
        )

        if self.compile_subproc.returncode != 0:
            raise ValueError('DART did not successfully compile.')


class DartCompile(object):
    """Class for a dart build = mkmf + compile and its resulting objects."""
    def __init__(
        self,
        source_dir: str,
        mkmf_template: str,
        work_dirs: list=['models/wrf_hydro/work'],
        mpi: bool=True, 
        build_dir: str = None,
        overwrite: bool = False
    ):

        self.source_dir = pathlib.PosixPath(source_dir).absolute()
        self.mkmf_template = self.source_dir / ('build_templates/' + mkmf_template)
        if type(work_dirs) is not list:
            work_dirs = [work_dirs]
        self.work_dirs = [pathlib.PosixPath(ww) for ww in work_dirs]
        self.mpi = mpi
        self.build_dir = build_dir
        self.overwrite = overwrite

        self.git_hash = get_git_revision_hash(self.source_dir)

        # mkmf establishment
        mkmf_dir = self.source_dir / 'build_templates'
        mkmf_target = mkmf_dir / 'mkmf.template'
        if mkmf_target.exists():
            mkmf_target.unlink()
        mkmf_target.symlink_to(self.mkmf_template)

        self.object_id = None
        """str: A unique id to join object to compile directory."""

        ## Setup directory paths
        if self.build_dir is not None:
            self.build_dir = pathlib.PosixPath(self.build_dir)
            self.build_dir.mkdir()
            # TODO(JLM): enforce that basename(build_dir) is experiment_dir

        for www in self.work_dirs:
            ww = str(www)
            dart_work = DartWork(
                ww,
                self.source_dir,
                self.build_dir,
                self.mpi
            )
            ww_repl = ww.replace('/','__')
            self.__dict__.update({ww_repl:dart_work})

        # Add in unique ID file to match this object to prevent assosciating
        # this directory with another object
        self.object_id = str(uuid.uuid4())
        with open(self.build_dir.joinpath('.uid'),'w') as f:
            f.write(self.object_id)

        self.pickle()
        print('DART successfully compiled into ' + str(self.build_dir))

    def pickle(self):
        with open(self.build_dir.joinpath('DartCompile.pkl'), 'wb') as f:
            pickle.dump(self, f, 2)


class HydroDartRun(object):
    """Class for dart and wrf-hydro runs (currently just filter?)."""
    def __init__(
        self,
        dart_compile: DartCompile,
        wrf_hydro_ens_sim, #: WrfHydroEnsembleRun,
        config: dict={}
    ):

        self.run_dir = pathlib.PosixPath(config['experiment']['run_dir'])
        """The absolute path to the hydro-dart run dir."""
        self.config = copy.deepcopy(config)
        """The configuation from the experiment setup."""

        self.wrf_hydro_ens_sim_pkl = self.run_dir / "WrfHydroEnsembleRun.pkl"
        self.exp_dir = self.run_dir / 'experiment_dir'
        self.dart_compile_pkl = self.exp_dir / (config['dart']['build_dir'] + "/DartCompile.pkl")

        self.binaries = {}
        self.scripts = {}

        # jobs_pending
        # job_active
        # jobs_completed
        self.pickle()        

    def advance_ensemble(
        self, 
        model_start_time: datetime=None,
        model_end_time: datetime=None,
        job_entry_cmd: str=None,
        job_exit_cmd: str=None,
        afterok: str=None,
        afterany: str=None,
        hold: bool=False
    ):

        # Setup job and scheduler.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Create a default job
            the_job = Job(
                nproc=self.config['run_experiment']['wrf_hydro_ens_advance']['nproc'],
                entry_cmd=job_entry_cmd,
                exit_cmd=job_exit_cmd
            )

            # TODO(JLM): add the entry and exit to the job

            # Create a dummy scheduler
            # This could be done with **{config['run_experiment']['wrf_hydro_ens_advance']['job_name']}
            if self.config['run_experiment']['wrf_hydro_ens_advance']['with_filter']:
                the_sched = None
                ppn_max = self.config['run_experiment']['dart']['scheduler']['ppn_max']
                nproc_model = \
                    self.config['run_experiment']['wrf_hydro_ens_advance']['nproc']
                # Leave a processor for the OS.
                n_mem_simul = math.floor((int(ppn_max) - 1) / int(nproc_model))
                print("n_mem_simul: ", n_mem_simul)
            else:
                the_sched = Scheduler(
                    job_name=self.config['run_experiment']['wrf_hydro_ens_advance']['job_name'],
                    account=self.config['run_experiment']['wrf_hydro_ens_advance']['account'],
                    nproc=self.config['run_experiment']['wrf_hydro_ens_advance']['nproc'],
                    nnodes=self.config['run_experiment']['wrf_hydro_ens_advance']['nnodes'],
                    walltime=self.config['run_experiment']['wrf_hydro_ens_advance']['walltime'],
                    queue=self.config['run_experiment']['wrf_hydro_ens_advance']['queue'],
                    email_who=self.config['run_experiment']['wrf_hydro_ens_advance']['email_who'],
                    email_when=self.config['run_experiment']['wrf_hydro_ens_advance']['email_when'],
                    afterok=afterok
                )

            # TODO(JLM): add afterany

        the_job.scheduler = the_sched

        ens_sim = pickle.load(open(self.wrf_hydro_ens_sim_pkl, 'rb'))
        # TODO (JLM): this is in place of the "sweeper" job or part of the submit script.
        ens_sim.collect_ensemble_runs(n_mem_simultaneous=n_mem_simul)

        if model_start_time is None:
            model_start_time = get_ens_dotfile_end_datetime(self.run_dir)
        if model_end_time is None:
            model_end_time = \
                model_start_time + \
                timedelta(hours=self.config['run_experiment']['time']['advance_model_hours'])

        the_job.model_start_time = model_start_time
        the_job.model_end_time = model_end_time

        ens_sim.add(the_job)
        ens_sim.run(
            hold=hold,
            n_mem_simultaneous=n_mem_simul
        )

        self.pickle()

    def pickle(
        self,
        path: pathlib.PosixPath=None
    ):
        if path is None:
            path = self.config['experiment']['experiment_dir']
        filepath = path / 'HydroDartRun.pkl' 

        with open(filepath, 'wb') as f:
            pickle.dump(self, f, 2)

    def establish_run_dir(self):

        dart_compile = pickle.load(open(self.dart_compile_pkl,'rb'))

        binary_src_path = dart_compile.build_dir
        the_binary = self.config['run_experiment']['dart']['exe']
        binary_src = binary_src_path / the_binary
        binary_link = self.config['experiment']['run_dir'] / the_binary
        binary_link.symlink_to(binary_src)
        self.binaries[the_binary] = binary_link

        for ss in dart_compile.run_scripts:
            script_src = ss
            script_base = os.path.basename(str(script_src))
            script_link = self.config['experiment']['run_dir'] / script_base
            script_link.symlink_to(script_src)
            self.scripts[script_base] = script_link
