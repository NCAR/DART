function [area] = heron(a, b, c)

s = (a + b + c) /2;

arg = (s * (s - a) * (s - b) * (s - c));

if(arg < -1e-10)
   fprintf('STOPPING IN HERON \n');
   [a b c s arg]
   stop
end
   
% Make sure we don't roundoff to a negative
if(arg <= 0) 
   area = 0;
   return
end

area = sqrt(arg);


