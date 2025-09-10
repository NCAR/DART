% Do some straightforward perturbing of aether ensemble files

% Plot values from aether neutrals block files on sphere
nblocks = 6;
ens_size = 10;

n_prior_base_name = 'B6_AETHER_INPUT_FILES/restartOut/neutrals_m';
n_post_base_name = 'B6_AETHER_OUTPUT_FILES/neutrals_m';
n_inc_base_name = './increments//neutrals_inc_m';
i_prior_base_name = 'B6_AETHER_INPUT_FILES/restartOut/ions_m';
i_post_base_name = 'B6_AETHER_OUTPUT_FILES/ions_m';
i_inc_base_name = './increments/ions_inc_m';

% Loop through the blocks and set values
for block = 0:nblocks -1 
   % Get the lat and lon info for this block
   block_prelim = int2str(10000 + block);
   block_final = block_prelim(2:5);

   for ens = 0:ens_size-1
      ens_prelim = int2str(10000 + ens);
      ens_final = ens_prelim(2:5);
      n_prior_file_name = strcat(n_prior_base_name, ens_final, '_g', block_final, '.nc');
      n_post_file_name = strcat(n_post_base_name, ens_final, '_g', block_final, '.nc');
      n_inc_file_name = strcat(n_inc_base_name, ens_final, '_g', block_final, '.nc');
      i_prior_file_name = strcat(i_prior_base_name, ens_final, '_g', block_final, '.nc');
      i_post_file_name = strcat(i_post_base_name, ens_final, '_g', block_final, '.nc');
      i_inc_file_name = strcat(i_inc_base_name, ens_final, '_g', block_final, '.nc');

      command_string = "rm " + n_inc_file_name;
      status = system(command_string);

      command_string = "ncdiff " + n_prior_file_name + " " + n_post_file_name + " " + n_inc_file_name;
      status = system(command_string);
      
      command_string = "rm " + i_inc_file_name;
      status = system(command_string);

      command_string = "ncdiff " + i_prior_file_name + " " + i_post_file_name + " " + i_inc_file_name;
      status = system(command_string);

   end
end

