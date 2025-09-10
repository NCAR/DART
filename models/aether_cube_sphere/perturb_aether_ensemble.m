% Do some straightforward perturbing of aether ensemble files

% Plot values from aether neutrals block files on sphere
nblocks = 6;
ens_size = 10;

g_base_name = 'B6_AETHER_INPUT_FILES/restartOut/grid_g';
n_base_name = 'B6_AETHER_INPUT_FILES/restartOut/neutrals_m';
i_base_name = 'B6_AETHER_INPUT_FILES/restartOut/ions_m';

% Loop through the blocks and set values
for block = 0:nblocks -1 
   % Get the lat and lon info for this block
   block_prelim = int2str(10000 + block);
   block_final = block_prelim(2:5);
   g_file_name = strcat(g_base_name, block_final, '.nc')

   % Get the lat and lon for more advanced perturbing
   % Read the grid file lat and lons
   lat = ncread(g_file_name, 'Latitude');
   lon = ncread(g_file_name, 'Longitude');

   % --------- Loop through the variables for the neutrals files --------
   % Get info about the variables from the first ensemble member 0000 file
   ens_prelim = int2str(10000 + 0);
   ens_final = ens_prelim(2:5);
   n_file_name_0 = strcat(n_base_name, ens_final, '_g', block_final, '.nc')
   i_file_name_0 = strcat(i_base_name, ens_final, '_g', block_final, '.nc')

   % Loop through all the variables in the neutrals files to perturb
   nc_info = ncinfo(n_file_name_0);
   % Number of total variables (including time)
   nvars = size(nc_info.Variables, 2);
  
   % Loop through the non-time variables 
   for ivar = 2:nvars
      var_name = nc_info.Variables(ivar).Name

      % Get the base variable from the first ensemble member, the current file name
      % Read variable for this variable from ensemble member 0
      var = ncread(n_file_name_0, var_name);
    
      % Loop through the rest of the ensemble members and perturb
      for ens = 1:ens_size - 1

         ens_prelim = int2str(10000 + ens);
         ens_final = ens_prelim(2:5);
         n_file_name = strcat(n_base_name, ens_final, '_g', block_final, '.nc')

         % Copy the ensemble member 0 file to the ensemble 1 using the shell; first ivar only
         if(ivar == 2) copyfile(n_file_name_0, n_file_name); end

         % Perturb the variable
         var_range = range(var, 'all')
         if(var_range == 0) 
            pert_var = var + 0.01 * var(1, 1, 1) * ens;
         else
            pert_var = var + var_range * 0.01 * ens;
         end
         % Write the variable
         ncwrite(n_file_name, var_name, pert_var);
      end
   end

   % --------- Loop through the variables for the ions files --------
   nc_info = ncinfo(i_file_name_0);
   % Number of total variables (including time)
   nvars = size(nc_info.Variables, 2);
  
   % Loop through the non-time variables 
   for ivar = 2:nvars
      var_name = nc_info.Variables(ivar).Name

      % Get the base variable from the first ensemble member, the current file name
      % Read variable for this variable from ensemble member 0
      var = ncread(i_file_name_0, var_name);
    
      % Loop through the rest of the ensemble members and perturb
      for ens = 1:ens_size - 1

         ens_prelim = int2str(10000 + ens);
         ens_final = ens_prelim(2:5);
         i_file_name = strcat(i_base_name, ens_final, '_g', block_final, '.nc')

         % Copy the ensemble member 0 file to the ensemble 1 using the shell
         if(ivar == 2) copyfile(i_file_name_0, i_file_name); end

         % Perturb the variable
         var_range = range(var, 'all')
         if(var_range == 0) 
            pert_var = var + 0.01 * var(1, 1, 1) * ens;
         else
            pert_var = var + var_range * 0.01 * ens;
         end
         % Write the variable
         ncwrite(i_file_name, var_name, pert_var);
      end
   end
end
