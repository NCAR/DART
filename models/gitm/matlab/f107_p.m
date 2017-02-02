function [f,d,fa]=f107_p(p,y,m,d1,nd,pl)

% F107_P( P, YYYY, MM, DD, ND ) - Plot the contents of f107.txt for a selected
% period of days, vs a minute-based x-axis.
%
% P is a string representing the path to f107.txt you want to plot
% YYYY is start year (as in 2002)
% MM is start month (as in 12, or if want to supply in number-of-days-since-beginning of year format, set to NAN)
% DD is start day (as in 15 or 335 (read all of previous line))
% ND is number of days you want to plot (as in 2)
%
% Example:
% plot(1:2880,150+rand(1,2880)) %fake data
% hold on
% f107_p('~/',2002,12,1,10)
% hold off

% DART $Id$
% CREDIT: Alexey Morozov

%% READING DATA PART
po=pwd;
fn ='f107.txt';
fn2='f107.txt.new';
cd(p)

[~,r] = unix(['wc -l ' fn]);
r=str2num(r(1:find(r>57,1)-1)); %take only the number part of reply
unix(['tail -n ' num2str(r-16) ' ' fn ' > ' fn2]); %remove text on the top
unix(['sed -i''.tmp'' ''s/"//g'' ' fn2]); %remove quotes
unix(['sed -i''.tmp'' ''s/-/ /g'' ' fn2]); %replace - with space
unix(['sed -i''.tmp'' ''s/:/ /g'' ' fn2]); %replace : with space

a=load(fn2);
cd(po)
yy=a(:,1);
mo=a(:,2);
dd=a(:,3);
hh=a(:,4); %f107.txt doesn't use anything but hh=0, so hh is ignored in the rest of this program
mm=a(:,5); % " "
f1=a(:,6);


dna=datenum(yy, mo, dd); %datenum available from f107.txt

if isnan(m) %user requested the number-of-days-since-beginning of year format
    dnr=datenum(y,0,0)+d1; %datenum requested
else
    dnr=datenum(y,m,d1); %datenum requested
end

f1a=0*f1;
n=40; %number of previous and next points to be averaged over (40 implies 81 day centered average)
for i=(n+1):(length(f1)-n)
    f1a(i) = mean(  f1( (i-n):(i+n) )  );
end

%% PLOTTING PART

i=find(dna==dnr);
f=f1(i:i+nd);
fa=f1a(i:i+nd);
d=(0:nd)*1440; %if you want linear minutes on x-axis
% d=(0:nd); %if you want linear days on x-axis

% d=d1c+(0:nd); %if you want linear days on x-axis
% ds=dd(i:i+nd);
if pl
    plot(d,f,'r')
    hold on
    plot(d,fa,'b')
    plot(d,(f+fa)/2,'m') %the actual value used in gitm -average of f107 and f107a
    hold off
    
    set(gca,'XTick',d) %if you selected linear days above, you can change XTickLabel to calendar days here
    % set(gca,'XTickLabel',ds)
    
    grid on
    box on
    legend('Daily F10.7','81-day central average of f107, called f107a','The value used in GITM, ie (f107+f107a)/2')
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
