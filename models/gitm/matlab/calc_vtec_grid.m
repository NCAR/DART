function VTecTS=calc_vtec_grid(AltT, RhoT, wlon, wlat)

% DART $Id$
% CREDIT: Alexey Morozov

nAlts=length(AltT);

IdseTS=squeeze(RhoT(wlon,wlat,:,:));
H=diag(ones(nAlts,1))+diag(ones(nAlts-1,1),1);
H=H(1:nAlts-1,:); %summation operator (needed for the trapezoidal rule - Integral = Sum_{i=1}^{N-1} (h_{i+1}-h_{i})/2*(f_{i+1}+f_{i}), so H*Idse gives you (f_{i+1}+f_{i})
IdseTSs=H*IdseTS; % (f_{i+1}+f_{i})
da=diff(AltT')/2; % (h_{i+1}-h_{i})/2
VTecTS=da*IdseTSs; % Integral = Sum_{i=1}^{N-1} (h_{i+1}-h_{i})/2 * (f_{i+1}+f_{i})
VTecTS=VTecTS*10^-16; %convert it into TECU

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
