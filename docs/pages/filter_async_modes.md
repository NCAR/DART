---
title: Filter Async
layout: default
---

# Filter async modes

## Options for parallelism both in DART and in the model advances:

<span id="async0"></span>  
Simplest case, async=0:  
![](../images/async0.gif)  
This is a single MPI executable, with each call to the model being
simply a subroutine call from each MPI task.  
  
[To the DART mpi intro document](dart_mpi.md#async0)  
  
<span id="async2"></span>  
Parallel advance, async=2:  
![](../images/async2a.gif)  
The filter executable is one MPI program, and the model is a single,
sequential executable. Each MPI task uses the unix "system()" call to
invoke a shell script (advance_model.csh) which runs the models as
independent programs.  
  
[To the DART mpi intro document](dart_mpi.md#async2)  
  
Other views of how the async=2 option is structured; these may be more
or less helpful.  
Parallel advance, async=2:  
![](../images/async2_v1.gif)  
  
Parallel advance, async=2, second version:  
![](../images/async2_v2.gif)  
  
Parallel model advance, async=2, showing how data is communicated
between filter and the model thru intermediate files. IC are 'initial
condition' files, UD are 'updated' files.  
![](../images/async2_wfiles.gif)  
    
<span id="async4"></span>  
Parallel model advance, async=4:  
![](../images/async4.gif)  
The filter executable is one MPI program, and the model is also an MPI
program. The filter executable communicates with the runme_filter shell
script, which sequentially invokes mpirun to advance each of the model
runs, one per ensemble member, still using advance_model.csh.  
  
[To the DART mpi intro document](dart_mpi.md#async4)  
   
Parallel model advance, async=4, showing how data is communicated
between filter and the model thru intermediate files. IC are 'initial
condition' files, UD are 'updated' files.  
![](../images/async4_wfiles.gif)  
  
-----
