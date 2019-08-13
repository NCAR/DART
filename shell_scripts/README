!---------------------------------------------
!> README file on how to use the bash scripts
!> scripts are written for lsf queuing system
!> they should be modified to work with others such as slurm

!> FESOM executables are called externally detached from DART
!> therefore no need for an advance model

environment.load           :: environment variables, relevant directories, experiment specification are set
                              here. There is a template called environment.template to be copied
                              here. This script is sourced by every other one. serial

ensemble.launch            :: main script modifying ens_members.template.lsf and calling ens_members.${EXPINFO}.lsf.
                              Keeps also an experiment specific brief information about the experiment which should
                              be modified before launching the scripts. serial

ens_members.template.lsf   :: calls initialize.template, forward_model.template check_ensemble.template subsequently.
                              serial job submitted to the queue

initialize.template        :: is called only once at the beginning of the experiment. Sets the experiment directory,
                              copies initial ensemble, namelists. serial job submitted to queue

forward_model.template     :: submits a job array including all ensemble members. parallel

check_ensemble.template    :: checks if the forwarding all members is finished. If so, first calls filter.template
                              and then calls finalize.template to conclude current assimilation cycle

filter.template            :: runs the filter after all members are forwarded to performs the analysis

finalize.template          :: checks if the whole experiment is finished. If so, stops. Otherwise, resubmits
                              ens_members.${EXPINFO}.lsf for the next assimilation cycle. for the next assimilation cycle.
