function time_plot(ttr,truth,latr, ...
    ttm,middle,latm, ...
    tms,meas,lams, ...
    tmd,mead,lamd, ...
    tmc,meac,meacsize,meaccol,lamc, ...
    ten,meanR,meanO,lamean, ...
    tes,sprR,sprO,laspr, ...
    teno,eno,laeno, ...
    distr,dsh,ddom,lad, ...
    leloc,leint,leord,...
    a, ...
    s )

% a custom plotter for the truth, measurement, ensemble, initial distribution (any of which can
% be ommitted via supplying an [] instead of values).
%
%ttr - Time for TRuth
%truth - is the true state
%latr - legend label for truth
%
%ttm - Time for Middle
%middle - is the middle state
%latm - legend label for middle
%
%tms - Time for MeaSurement
%meas - MEASurement (ie usually = truth+noise)
%lams - legend label for meas
%
%tmd - Time for Measurement spread (+-standard Deviation)
%mead - MEAsurement standard Deviation
%lamd - legend label for mead
%
%tmc - Time for Measurement that is sCattered
%meac - MEAsurement sCattered (ie satellite flown through truth)
%meacsize - size of scatter circles (50 is ok)
%meaccolor - length(tms) by 3 matrix specifying colors of each scatter circle
%lamc - legend label for meac
%
%ten - Time for ENsemble means and spreads
%meanR - pRior ensemble mean
%meanO - pOsterior ens. mean
%lamean - legend label for mean
%
%tes - Time for Ensemble Spread (should be = to ten, but if [], then spread is not plotted)
%sprR - pRior spread (+-standard deviation over the prior ensemble)
%sprO - pOsterior spread (+-standard deviation over the posterior ensemble)
%laspr - legend label for spread
%
%teno - Time for Ensemble Number of Observations
%eno - number of observations used
%laeno - legend label
%
%distr - initial state DISTRibution (ie [mean, var] = 2 element vector)
%dsh - how to scale and shift the distr so it would look nice (not to scale)
%ddom - the domain over which to draw the distribution (ie -3:.1:3)
%lad - legend label for distr
%
%leloc - LEgend LOCation ('NE' or 'Best' or ...) if is set to NAN, then no legend will be displayed whatsoever
%leint - LEgend INTerpreter ('tex' (normal) or 'latex' (fancy))
%leord - LEgend ORDer: standard is 'meas truth mean spr meac dist', which is equivalent to leord=[1 2 3 4 5 6], which is defaulted to if leord = []. More reasonable might be [2 1 3 4 5 6] etc (if you enabled all of them, if not, you have to exclude some and keep track of the order).
%
%a - alpha - transparency of patches, (0 to 1), default is 1
%
%s - scaling factor. If you want to deal with Rho, you might wanna set it
%to 10^12, default is 1

% DART $Id$
% CREDIT: Alexey Morozov

%% truth
%no preprocessing to be done here

%% measurement
%no preprocessing to be done here

%% ensemble
if ~isempty(ten)
    tenRO  = nan(1,2*length(ten)); %time for the sawtooth (need to duplicate it)
    meanRO = nan(1,2*length(ten)); %mean-sawtooth (need to mesh prior and posterior together: [R1 O1 R2 O2 R3 O3 ...]
    
    tenRO(1:2:end)=ten(:)';
    tenRO(2:2:end)=ten(:)';
    meanRO(1:2:end)=meanR(:)';
    meanRO(2:2:end)=meanO(:)';
end

if ~isempty(tes)
    sprRO  = nan(1,2*length(ten)); %spr-sawtooth (see comments above)
    sprRO(1:2:end)=sprR(:)';
    sprRO(2:2:end)=sprO(:)';
end
%% distribution
if ~isempty(distr)
    di = 1 / sqrt( 2*pi*distr(2) ) * exp( -(ddom-distr(1)).^2 / (2*distr(2)) ) ;
end

%% plotting
cla
set(gca,'position',[.13 .11 .78 0.82])
hold on

h=[]; %array of handles to plots
la={}; %cell array of legend entries for the plots in h

%distribution
if ~isempty(distr)
    hdi = patch( di/max(di)*dsh(1)+dsh(2), s*ddom, 'g','edgecolor',[0 .5 0],'FaceAlpha',a); %,'FaceAlpha',0.5
    % di/max(di)*dsh(1)+dsh(2) makes it not to scale, but helps to make it more visible
    h=[h hdi];
    la{1,end+1}=lad;
end

%champ localized dots
if ~isempty(tmc)
    hmc = scatter(tmc,s*meac,meacsize,meaccol,'filled');
    h=[h hmc];
    la{1,end+1}=lamc;
end

%ensemble number of obs used
if ~isempty(teno)
    heno=plot(teno,eno,'b+');
    h=[h heno];
    la{1,end+1}=laeno;
end

%measurement spread
if ~isempty(tmd)
    hmd=patch([tmd(:)' fliplr(tmd(:)')],s*[meas(:)'+mead(:)' fliplr(meas(:)'-mead(:)')],'r','facecolor',[255 153 153]/255,'edgecolor','r','FaceAlpha',a,'linestyle','--');%,'FaceAlpha',0.5
    h=[h hmd];
    la{1,end+1}=lamd;
end

%ensemble spread
if ~isempty(tes)
    hspr=patch([tenRO(:)' fliplr(tenRO(:)')],s*[meanRO(:)'+sprRO(:)' fliplr(meanRO(:)'-sprRO(:)')],'b','facecolor',[135 206 250]/255,'edgecolor','b','FaceAlpha',a);%,'FaceAlpha',0.5
    h=[h hspr];
    la{1,end+1}=laspr;
end




%ensemble mean
if ~isempty(ten)
    hmean=plot(tenRO,s*meanRO,'b','linewidth',2);
    h=[h hmean];
    la{1,end+1}=lamean;
end



%measurement
if ~isempty(tms)
    hms=plot(tms,s*meas,'--r','linewidth',2);
    h(1,end+1)=hms;
    la{1,end+1}=lams;
end


%middle
if ~isempty(ttm)
    htm=plot(ttm,s*middle,'-.k','linewidth',1);
    h(1,end+1)=htm;
    la{1,end+1}=latm;
end

%truth
if ~isempty(ttr)
    htr=plot(ttr,s*truth,':k','linewidth',2);
    h(1,end+1)=htr;
    la{1,end+1}=latr;
end



box on
grid on
hold off

%% legend
if ~isnan(leloc)
    h=fliplr(h); %want truth displayed first in legend and distr last in legend (but plotted in the reverse order, so that truth is front-most and distr is back-most)
    la=fliplr(la); %want -//-
    
    if ~isempty(leord)
        h=h(1,leord); %if user supplied a fancy order, legend will have that order
        for i=1:length(leord)
            la1{1,i}=la{1,leord(i)};
        end
        la=la1;
    end
    
    le=legend(h,la,'location',leloc);
    set(le,'interpreter',leint)
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
