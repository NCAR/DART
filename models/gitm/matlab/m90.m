function [to,yo]=m90(ti,yi,m)

%calculate 90 minute averages
%ti - time in (hours?)
%yi - values in
%m - number of (units of ti) to calculate averages over
%to - time out
%yo - values out

% DART $Id$
% CREDIT: Alexey Morozov

to=ti(1):m:ti(end); 
yo=0*to; %preallocate

for i=2:length(to)
    cur=find( to(i-1)<ti & ti<=to(i));
    yo(i)=sum(yi(cur))/length(cur);
end
yo(1)=yo(2);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
