function time_plot( ...
    ttr,truth,latr, ...
    tms,meas,lams, ...
    tmc,meac,meacsize,meaccol,lamc, ...
    ten,meanR,meanO,lamean, ...
    tes,sprR,sprO,laspr, ...
    distr,dsh,ddom,lad, ...
    leloc,leint,leord)

% a custom plotter for the truth, measurement, ensemble, initial distribution (any of which can
% be ommitted via supplying an [] instead of values).
%
%ttr - Time for TRuth
%truth - is the true state
%latr - legend label for truth
%
%tms - Time for MeaSurement
%meas - MEASurement (ie usually = truth+noise)
%lams - legend label for meas
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
%sprR - pRior spread (standard deviation over the ensemble)
%sprO - pOsterior spread (standard deviation over the ensemble)
%laspr - legend label for spread
%
%distr - initial state DISTRibution (ie [mean, var] = 2 element vector)
%dsh - how to scale and shift the distr so it would look nice (not to scale)
%ddom - the domain over which to draw the distribution (ie -3:.1:3)
%lad - legend label for distr
%
%leloc - LEgend LOCation ('NE' or 'Best' or ...) if is set to NAN, then no legend will be displayed whatsoever
%leint - LEgend INTerpreter ('tex' (normal) or 'latex' (fancy))
%leord - LEgend ORDer: standard is 'meas truth mean spr meac dist', which is equivalent to leord=[1 2 3 4 5 6], which is defaulted to if leord = []. More reasonable might be [2 1 3 4 5 6] etc (if you enabled all of them, if not, you have to exclude some and keep track of the order).

% DART $Id$
% CREDIT: Alexey Morozov

%% truth
%no preprocessing to be done here

%% measurement
%no preprocessing to be done here

%% ensemble
tenRO  = nan(1,2*length(ten)); %time for the sawtooth (need to duplicate it)
meanRO = nan(1,2*length(ten)); %mean-sawtooth (need to mesh prior and posterior together: [R1 O1 R2 O2 R3 O3 ...]
sprRO  = nan(1,2*length(ten)); %spr-sawtooth (see line above)
tenRO(1:2:end)=ten;
tenRO(2:2:end)=ten;
meanRO(1:2:end)=meanR;
meanRO(2:2:end)=meanO;
sprRO(1:2:end)=sprR;
sprRO(2:2:end)=sprO;

%% distribution
if ~isempty(distr)
    di = 1 / sqrt( 2*pi*distr(2) ) * exp( -(ddom-distr(1)).^2 / (2*distr(2)) ) ;
end

%% plotting
cla
hold on

h=[]; %array of handles to plots
la={}; %cell array of legend entries for the plots in h


if ~isempty(distr)
    hdi = patch( di/max(di)*dsh(1)+dsh(2), ddom, .55*[1,1,1],'edgecolor',.45*[1,1,1],'linestyle','--'); % di/max(di)*dsh(1)+dsh(2) makes it not to scale, but helps to make it more visible
    h=[h hdi];
    la{1,end+1}=lad;
end

if ~isempty(tmc)
    hmc = scatter(tmc,meac,meacsize,meaccol,'filled');
    h=[h hmc];
    la{1,end+1}=lamc;
end

if ~isempty(tes)
    hspr=patch([tenRO fliplr(tenRO)],[meanRO+sprRO fliplr(meanRO-sprRO)],.8*[1,1,1],'edgecolor',.5*[1,1,1]);
    h=[h hspr];
    la{1,end+1}=laspr;
end

if ~isempty(ten)
    hmean=plot(tenRO,meanRO,'color',.5*[1,1,1],'linewidth',2);
    h=[h hmean];
    la{1,end+1}=lamean;
end

if ~isempty(ttr)
    htr=stairs(ttr,truth,':k','linewidth',2);
    h(1,end+1)=htr;
    la{1,end+1}=latr;
end

if ~isempty(tms)
    hms=stairs(tms,meas,'--','color',.2*[1,1,1],'linewidth',2);
    h(1,end+1)=hms;
    la{1,end+1}=lams;
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

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
