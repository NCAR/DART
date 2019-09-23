
![DARTlogo](docs/images/Dartboard7.png)

## Welcome to DART, the Data Assimilation Research Testbed.

[quick-start instructions for the impatient](#QuickStart)

Extensive on-line documentation is available at the GitHub project web pages:
[https://ncar.github.io/DART](https://ncar.github.io/DART)
or in the repository at [docs/index.md](docs/index.md). It's probably a good idea
to be familiar with the [GitHub workflow](https://guides.github.com/introduction/flow/).

<table><tr>
<td><img height=250, src="./docs/images/DARTspaghettiSquare.gif"></td> 
<td><img height=250, src="./docs/images/science_nuggets/AssimAnim.gif"></td>
</tr></table>

[Extensive local documentation is included with DART.](docs/index.md)

A Matlab-based introduction is in the `docs/DART_LAB` directory.  
There are a set of PDF presentations along with hands-on Matlab exercises.  
This starts with a very basic introduction to data assimilation and covers 
several fundamental algorithms in the system.
A more exhaustive tutorial for data assimilation with DART is in PDF format at 
[docs/tutorial](docs/tutorial/index.html).

The DART Manhattan release documentation is on the web:
http://www.image.ucar.edu/DAReS/DART/Manhattan/documentation/html/Manhattan_release.html
and in the repository at:
[docs/html/Manhattan_release.html](docs/html/Manhattan_release.html)

There is a mailing list where we summarize updates to the DART repository
and notify users about recent bug fixes.
It is not generally used for discussion so it's a low-traffic list.
To subscribe to the list, click on
[Dart-users](http://mailman.ucar.edu/mailman/listinfo/dart-users).
If you use WRF, there is also a
[Wrfdart-users](http://mailman.ucar.edu/mailman/listinfo/wrfdart-users).

The Manhattan release is new and currently supports only a subset of the 
models.  We will port over any requested model so contact us if yours
is not on the list.  In the meantime, we suggest you check out our
'classic' release of DART which is the Lanai release plus additional
development features.  All new development will be based on the
Manhattan release but the 'classic' release will remain for those
models which already have the necessary features.

Contact us for more help or for more information on other models already
using DART or for how to add your model or observation types.

Thank you -
The DART Development Team.
dart at ucar.edu

## DART source code tree

The top level DART source code tree contains the following directories and files:

| Directory            | Purpose  |
| :--------------      | :------- |
| `assimilation_code/` | Low-level library and fonts required by NCAR Graphics and NCL |
| `build_templates/`   | Configuration files for installation |
| `developer_tests/`   | regression testing |
| `diagnostics/`       | routines to diagnose assimilation performance |
| `docs/`              | General documentation and DART_LAB tutorials |
| `models/`            | the interface routines for the models |
| `observations/`      | routines for converting observations and forward operators |
| **Files**            | **Purpose** |
| `CHANGELOG.md`       | Brief summary of recent changes |
| `COPYRIGHT.md`       | terms of use and copyright information |
| `README.md`          | Basic Information about DART |

## Bug reports and feature requests

Use the GitHub [issue tracker](https://github.com/NCAR/DART/issues) 
to submit a bug or request a feature.

## Citing DART

Cite DART using the following text:

> The Data Assimilation Research Testbed (Version X.Y.Z) [Software]. (2019). Boulder, Colorado: UCAR/NCAR/CISL/DAReS.  http://doi.org/10.5065/D6WQ0202

Update the DART version and year as appropriate.

---

<a name="QuickStart"></a>
## Quick-start for the impatient:

There are several large files that are needed to run some of the tests and examples but are not included
in order to keep the repository as small as possible. If you are interested in running *bgrid_solo*, 
*cam-fv*, or testing the *NCEP/prep_bufr* observation converter, you will need these files.
These files are available at:

| Release                | Size   | Filename  |
| :--------------        | :----- | :-------- |
| "Manhattan"        |  189M  | [Manhattan_large_files.tar.gz](https://www.image.ucar.edu/pub/DART/Release_datasets/Manhattan_large_files.tar.gz) |
| "wrf-chem.r13172"  |  141M  | [wrf-chem.r13172_large_files.tar.gz](https://www.image.ucar.edu/pub/DART/Release_datasets/wrf-chem.r13172_large_files.tar.gz) |
| "Lanai"            |  158M  | [Lanai_large_files.tar.gz](https://www.image.ucar.edu/pub/DART/Release_datasets/Lanai_large_files.tar.gz) |
| "Kodiak"           |  158M  | [Kodiak_large_files.tar.gz](https://www.image.ucar.edu/pub/DART/Release_datasets/Kodiak_large_files.tar.gz) |
| "Jamaica"          |   32M  | [Jamaica_large_files.tar.gz](https://www.image.ucar.edu/pub/DART/Release_datasets/Jamaica_large_files.tar.gz) |
| "Hawaii"           |   32M  | [Hawaii_large_files.tar.gz](https://www.image.ucar.edu/pub/DART/Release_datasets/Hawaii_large_files.tar.gz) |

Download the appropriate tar file and untar it into your DART repository. Ignore any warnings about
```tar: Ignoring unknown extended header keyword``` .

Go into the `build_templates` directory and copy over the closest
`mkmf.template`._compiler.system_ file into `mkmf.template`.

Edit it to set the NETCDF directory location if not in `/usr/local`
or comment it out and set $NETCDF in your environment.  *This NetCDF 
library must have been compiled with the same compiler
that you use to compile DART and must include the F90 interfaces.*

Go into `models/lorenz_63/work` and run *quickbuild.csh*.

> cd models/lorenz_63/work  
> ./quickbuild.csh  

If it compiles, *:tada:*!  Run this series of commands to do a very basic test:

> ./perfect_model_obs  
> ./filter  

If that runs, *:tada:* again!  Finally, if you have Matlab installed on
your system add '$DART/diagnostics/matlab' to your matlab search path 
and run the 'plot_total_err' diagnostic script while in the 
`models/lorenz_63/work` directory.  If the output plots and looks 
reasonable (error level stays around 2 and doesn't grow unbounded) 
you're great!  Congrats.

If you are planning to run one of the larger models and want to
use the Lorenz 63 model as a test, run ```./quickbuild.csh -mpi```.
It will build filter and any other MPI-capable executables with MPI.
*The 'mpif90' command you use must have been built with the same 
version of the compiler as you are using.*

If any of these steps fail or you don't know how to do them, go to the
DART project web page listed above for very detailed instructions that
should get you over any bumps in the process.

---
