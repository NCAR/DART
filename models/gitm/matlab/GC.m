function y=GC(c,z)

% DART $Id$
% CREDIT: Alexey Morozov

cz1=[-1/4 1/2 5/8 -5/3 0 1]; %in powers of (abs(z)/c), abs(z)>=0 & abs(z)<=c 
cz2=[1/12 -1/2 5/8 5/3 -5 4 -2/3]; %/(abs(z)/c) and abs(z)>=c & abs(z)<=2c
cz3=0; %abs(z)>=2c

y=zeros(size(z));

test=find(abs(z)>=0 & abs(z)<=c);
z1=z(test);
y(test) = polyval(cz1,abs(z1)/c);

test=find(abs(z)>=c & abs(z)<=2*c);
z1=z(test);
y(test) = polyval(cz2,abs(z1)/c)./(abs(z1)/c);

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
