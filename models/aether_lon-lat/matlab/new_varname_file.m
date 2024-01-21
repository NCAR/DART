function new_varname_file(data_dir, member, nblocks)

% Converts Aether restart file names to updated (2024-1-17) versions.
% Copy the grid files, for completeness, into a new directory.
% Run in that directory.
% Gets all the contents of the existing files
% and writes them to a new file, with new variable names.
% > new_varname_file(data_dir, member, nblocks)
%
% Files are version=2,netcdf=4.9.0,hdf5=1.12.2

% DAI/Aether/Aaron_names/restartOut.Sphere.1member
% neut_old = ["N_4S" "O2" "N2" "NO" "He" "N_2D" "N_2P" "H" "O_1D" "CO2" ...
% neut_new = ["N"    "O2" "N2" "NO" "He" "N_2D" "N_2P" "H" "O_1D" "CO2" ...
% testdata2 (ensemble from Aaron)
neut_old = ["O" "N" "O2" "N2" "NO" "He" "N_2D" "N_2P" "H" "O_1D" "CO2" ...
            "Temperature" "Zonal Wind"    "Meridional Wind" "Vertical Wind"]
neut_new = ["O" "N" "O2" "N2" "NO" "He" "N_2D" "N_2P" "H" "O_1D" "CO2" ...
            "Temperature" "velocity_east" "velocity_north"  "velocity_up"]
% neut_old = ["N_4S" "Zonal Wind"];
% neut_new = ["N"    "velocity_east"];

% ions_old = ["O+2P"  "Temperature (bulk ion)" ];
% ions_new = ["O+_2P" "Temperature_bulk_ion"   ];
% DAI/Aether/Aaron_names/restartOut.Sphere.1member
% ions_old = ["O+" "O2+" "N2+" "NO+" "N+" "He+" "O+2D"  "O+2P"  ...
% ions_new = ["O+" "O2+" "N2+" "NO+" "N+" "He+" "O+_2D" "O+_2P" ...
% testdata2
ions_old = ["O+" "O2+" "N2+" "NO+" "N+" "He+" "O+_2D" "O+_2P"  ...
            "Temperature (bulk ion)" "Temperature (electron)"]
ions_new = ["O+" "O2+" "N2+" "NO+" "N+" "He+" "O+_2D" "O+_2P" ...
            "Temperature_bulk_ion"   "Temperature_electron"]

global fname_old  fname
format compact

for b = 0:nblocks-1
   % Neutrals
   fname = sprintf('neutrals_m%04d_g%04d.nc', member, b)
   fname_old = strcat(data_dir,fname)

   create_file_skel()
   add_vars(neut_old, neut_new, b)

   % Ions
   fname = sprintf('ions_m%04d_g%04d.nc', member, b)
   fname_old = strcat(data_dir,fname)

   create_file_skel()
   add_vars(ions_old, ions_new, b)

end

%- - - - 
function create_file_skel()

% Create file and define dimensions

   global fname_old  fname ncid_old ncid_new

   ncid_old = netcdf.open(fname_old,'NOWRITE');
%    ncdisp(fname_old)
   [ndims_old,nvars_old,ngatts_old,unlimdimid] = netcdf.inq(ncid_old);

% I can't make CLOBBER work with whatever permissions Mac and umask allow,
% so remove any existing file manually.  Pathetic.
   system(['rm ',fname]);
% Ah, I wish they'd mentioned this in the .create page:
% Add write permission to this directory.
% NOPE, needs to operate on an existing file, which I can't create.  Brilliant.
%    fileattrib('.',"+w")
   cmode = netcdf.getConstant('NETCDF4');
   cmode = bitor(cmode,netcdf.getConstant('CLOBBER'));
   ncid_new = netcdf.create(fname,cmode);

% Get dimensions and write them to the new file.
   for d = 1:ndims_old
      [dimname, dimlen] = netcdf.inqDim(ncid_old,d-1);
      dimid = netcdf.inqDimID(ncid_old,dimname);
      dimid = netcdf.defDim(ncid_new,dimname,dimlen);
   end
   netcdf.endDef(ncid_new);

%- - - - 
function add_vars(vars_old, vars_new, b)

   global fname_old  fname ncid_old ncid_new

% time differs from all the others
   data = ncread(fname_old,"time");
   dim_list = {"time"};
   nccreate(fname,"time", Dimensions=dim_list, Datatype="double")
   ncwrite(fname,"time",data)

   for n = 1:length(vars_old)
      add_var(vars_old(n),vars_new(n), b)

      % Ions; add associated variables
      if  contains(fname,"ions") & ...
         ~contains(vars_old(n),"bulk") & ...
         ~contains(vars_old(n),"electron") 
         add_assoc_vars(vars_old(n), vars_new(n), b)
      end
   end
%    
   netcdf.close(ncid_new)
   netcdf.close(ncid_old)

% - - - -
function add_var(var_old, var_new, b)

   global fname_old  fname 

   if b == 0
      sprintf('Renaming %s to %s',var_old,var_new)
   end

   data = ncread(fname_old,var_old);
   att  = ncreadatt(fname_old,var_old,"units");

   dim_list = {"z","y","x"};
   nccreate  (fname, var_new, Dimensions=dim_list, Datatype="single")
   ncwrite   (fname, var_new, data)
   ncwriteatt(fname, var_new,"units",att)

%- - - - 
function add_assoc_vars(var_old, var_new, b)

% Variables with names associated with ions.
% example   'Parallel Ion Velocity (Zonal) (O+2P)'      ...
% The 'Temperature' part of the names is the same, but other parts are different,
% NOTE: These names have the \s removed, but Matlab+NetCDF puts them back in
%       in the new file.
   i_assoc_old = [ ...
      "Temperature"                        ...
      "Parallel Ion Velocity (Zonal)"      ...
      "Parallel Ion Velocity (Meridional)" ...
      "Parallel Ion Velocity (Vertical)"   ...
      "Perp. Ion Velocity (Zonal)"         ...
      "Perp. Ion Velocity (Meridional)"    ...
      "Perp. Ion Velocity (Vertical)"      ];
   i_assoc_new = [              ...
      "Temperature"             ...
      "velocity_parallel_east"  ...
      "velocity_parallel_north" ...
      "velocity_parallel_up"    ...
      "velocity_perp_east"      ...
      "velocity_perp_north"     ...
      "velocity_perp_up"        ];
 
   for a = 1:7
      i_old = strcat(i_assoc_old(a),' (',var_old,')');
      i_new = strcat(i_assoc_new(a),' (',var_new,')');
      add_var(i_old, i_new, b)
   end
