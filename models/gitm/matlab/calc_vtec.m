function VTecTS=calc_vtec(tt, LonT, LatT, AltT, RhoT, ...
    LonTS, LatTS)

% DART $Id$
% CREDIT: Alexey Morozov

nAlts=length(AltT);

IdseTS=interpn(LonT,LatT,AltT,tt,RhoT(:,:,:,:),...
    repmat(LonTS,nAlts,1),repmat(LatTS,nAlts,1),repmat(AltT,1,length(tt)),repmat(tt,nAlts,1), ...
    'linear');
H=diag(ones(nAlts,1))+diag(ones(nAlts-1,1),1);
H=H(1:nAlts-1,:); %summation operator (needed for the trapezoidal rule - Integral = Sum_{i=1}^{N-1} (h_{i+1}-h_{i})/2*(f_{i+1}+f_{i}), so H*Idse gives you (f_{i+1}+f_{i})
IdseTSs=H*IdseTS; % (f_{i+1}+f_{i})
da=diff(AltT')/2; % (h_{i+1}-h_{i})/2
VTecTS=da*IdseTSs; % Integral = Sum_{i=1}^{N-1} (h_{i+1}-h_{i})/2 * (f_{i+1}+f_{i})
VTecTS=VTecTS*10^-16; %convert it into TECU

if sum(isnan(VTecTS))/length(VTecTS)>0.1 
    error('more than 10% of the data is NaNs')
end

for i=1:length(VTecTS)
    if isnan(VTecTS(i))
        %     try
        j=find(~isnan(VTecTS(i:end)),1)+i-1;%find the next non-nan number
        VTecTS(i)=interp1([tt(i-1) tt(j)], [VTecTS(i-1) VTecTS(j)], tt(i));
        %     catch
        %must have started with a nan or the sequence ends with a nan. not
        %sure what to do. extrap?
        %     end
    end
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
