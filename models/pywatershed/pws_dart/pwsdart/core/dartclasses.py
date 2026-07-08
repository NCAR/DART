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

        static_file = path_dart / self.path_rel / 'static_file_list.txt'
        if static_file.exists():
            shutil.copy(static_file, build_dir / self.path_rel / 'static_file_list.txt')
        

        self.compile(path_dart, path_rel, mpi)

        # list the mkmfs and get the binaries.
        #mkmfs = list((pathlib.PosixPath(path_dart) / path_rel).glob('mkmf_*'))
        #dart_exes = [pathlib.PosixPath(str(mm).replace('/mkmf_','/')) for mm in mkmfs]
        dart_exe_names = [ "closest_member_tool", "filter", "model_mod_check", "perfect_model_obs",
                           "advance_time", "create_fixed_network_seq", "create_obs_sequence",
                           "fill_inflation_restart", "obs_common_subset", "streamflow_obs_diag",
                           "obs_sequence_tool", "obs_seq_to_netcdf", "obs_diag",
                           "create_identity_streamflow_obs", "convert_streamflow"
                         ]

        dart_exes = []
        for exe_name in dart_exe_names:
            exe_path = pathlib.Path(path_dart) / path_rel / exe_name
            if exe_path.exists():
                dart_exes.append(exe_path)
            else:
                print(f"[WARNING] Missing expected executable: {exe_path}")

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
        if mpi:
            build_cmd = './quickbuild.sh'
        else:
            build_cmd = './quickbuild nompi'

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
        work_dirs: list=['models/pywatershed/work'],
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

        #self.git_hash = get_git_revision_hash(self.source_dir)

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


