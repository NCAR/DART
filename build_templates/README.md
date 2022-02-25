
#### Building and Customizing the 'mkmf.template' file

A series of templates for different compilers/architectures can be found in
the *DART/build_templates* directory and have names with extensions
that identify the compiler, the architecture, or both. This is how you
inform the build process of the specifics of your system. **Our intent
is that you copy one that is similar to your system into
`DART/build_templates/mkmf.template` and customize it.**

Go into the `build_templates` subdirectory and copy over the closest `mkmf.template.<compiler system>` file into `mkmf.template`.

Edit `mkmf.template` to set the `NETCDF` directory location if not in
`/usr/local` or comment it out and set `$NETCDF` in your environment.

For more information on how to customize the machine-specific resources,
[see here](https://dart.ucar.edu/pages/Getting_Started.html#customizations).
