src_dir: ../.
exclude_dir: ../docs/_api
             ../docs/_site/api
             ../models/bgrid_solo/fms_src
             ../observations/obs_converters/NCEP
             ../assimilation_code/programs/system_simulation
exclude: mpisetup.f90
         obs_seq.F
output_dir: ./ford_output
page_dir: ./pages
project: DART
project_github: https://github.com/NCAR/DART
project_website: https://ncar.github.io/DART
summary: **D**ata **A**ssimilation **R**esearch **T**estbed
author: DAReS
author_description: The NCAR Data Assimilation Research Section
display: public
source: false
graph: true
graph_maxdepth: 3
coloured_edges: true
search: true
md_extensions: markdown.extensions.toc
preprocessor: gfortran
css: ./_api/ncar.css

Click Documentation for release notes, getting started, and other documentation
