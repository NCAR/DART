% Do some straightforward perturbing of aether ensemble files
% found in the TEST_INPUT directory
% The perturbations are all constant scaled fields; resulting correlations are all 1.

% Number of blocks for default aether restart files
nblocks = 6;
% Size of ensemble to create
ens_size = 10;

% Directory containing files to be perturbed
pert_dir = 'TEST_INPUT';
g_base_name = strcat(pert_dir, '/grid_g');
n_base_name = strcat(pert_dir, '/neutrals_m');
i_base_name = strcat(pert_dir, '/ions_m');

% Loop through the blocks
for block = 0:nblocks -1 
   block_prelim = int2str(10000 + block);
   block_final = block_prelim(2:5);
   g_file_name = strcat(g_base_name, block_final, '.nc');

   % Get the lat and lon for more advanced perturbing
   lat = ncread(g_file_name, 'Latitude');
   lon = ncread(g_file_name, 'Longitude');

   ens_prelim = int2str(10000 + 0);
   ens_final = ens_prelim(2:5);
   % File names for neutral and ions ensemble member 0
   n_file_name_0 = strcat(n_base_name, ens_final, '_g', block_final, '.nc');
   i_file_name_0 = strcat(i_base_name, ens_final, '_g', block_final, '.nc');

   % Loop through all the variables in the neutrals files to perturb
   nc_info = ncinfo(n_file_name_0);
   % Number of total variables (including time)
   nvars = size(nc_info.Variables, 2);
  
   % Loop through the non-time variables 
   for ivar = 2:nvars
      var_name = nc_info.Variables(ivar).Name;

      % Get the base variable from the first ensemble member, the current file name
      var = ncread(n_file_name_0, var_name);
    
      % Loop through the rest of the ensemble members and perturb
      for ens = 1:ens_size - 1

         ens_prelim = int2str(10000 + ens);
         ens_final = ens_prelim(2:5);
         n_file_name = strcat(n_base_name, ens_final, '_g', block_final, '.nc');

         % Copy the ensemble member 0 file to the ensemble ens using the shell; first ivar only
         if(ivar == 2) copyfile(n_file_name_0, n_file_name); end

         % Compute a normalized range for perturbation size
         if(block == 0) 
            neutrals_var_range(ivar) = range(var, 'all');
            if(neutrals_var_range(ivar) == 0) 
               neutrals_var_range(ivar) = var(1, 1, 1);
            end
         end

         % Perturb the variable
         pert_var = var + neutrals_var_range(ivar) * 0.01 * ens;

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
      var_name = nc_info.Variables(ivar).Name;

      % Get the base variable from the first ensemble member, the current file name
      var = ncread(i_file_name_0, var_name);
    
      % Loop through the rest of the ensemble members and perturb
      for ens = 1:ens_size - 1

         ens_prelim = int2str(10000 + ens);
         ens_final = ens_prelim(2:5);
         i_file_name = strcat(i_base_name, ens_final, '_g', block_final, '.nc');

         % Copy the ensemble member 0 file to the ensemble 1 using the shell
         if(ivar == 2) copyfile(i_file_name_0, i_file_name); end

         % Compute a normalized range for perturbation size
         if(block == 0) 
            ions_var_range(ivar) = range(var, 'all');
            if(ions_var_range(ivar) == 0) 
               ions_var_range(ivar) = var(1, 1, 1);
            end
         end

         % Perturb the variable
         pert_var = var + ions_var_range(ivar) * 0.01 * ens;

         % Write the variable
         ncwrite(i_file_name, var_name, pert_var);
      end
   end
end
