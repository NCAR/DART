%%
%

file2 = '/glade/scratch/thoar/tiegcm_test/test2/Posterior_Diag.nc';
file4 = '/glade/scratch/chihting/job_test/test4/Posterior_Diag.nc';

vars = {'NE','TN','O1','O2','UN','VN','OP'};

for ivar = 1:length(vars)

   v2 = nc_varget(file2, vars{ivar});
   v4 = nc_varget(file4, vars{ivar});
   
   copy = 1; % Just checking the ensemble mean
   
   for iz = 1:size(v2,2)
      org = squeeze(v2(copy,iz,:,:));
      new = squeeze(v4(copy,iz,:,:));
      change = new-org;
  
      dmin = min([org(:); new(:)]);
      dmax = max([org(:); new(:)]);
 
      subplot(3,1,1);
         imagesc(org, [dmin dmax]); set(gca,'YDir','normal')
         title(sprintf('Original %s at level %d',vars{ivar}, iz))
         colorbar
   
      subplot(3,1,2);
         imagesc(new); set(gca,'YDir','normal')
         title(sprintf('Proposed %s at level %d',vars{ivar}, iz))
         colorbar
         
      subplot(3,1,3);
         imagesc(change); set(gca,'YDir','normal')
         title(sprintf('Difference (new-old) at level %d',iz))
         colorbar

      % Find the relative change in orders of magnitude and print if any are 'big'

      maxdiff   =  max(abs(change(:)));
      indx      = find(abs(change(:)) == maxdiff);
      relchange = log10(org(indx)) - log10(maxdiff);

      if ( relchange < 5 )
         fprintf('%s level %d has some non-neglible changes ... %f\n',vars{ivar}, iz, relchange)
      end
   
      fprintf('Showing %s level %d of %d ... pausing\n',vars{ivar}, iz, size(v2,2))
      pause(1.0)  % in seconds.
   end
end 
