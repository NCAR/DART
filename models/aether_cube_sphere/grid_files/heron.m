function [area] = heron(a, b, c)

s = (a + b + c) /2;

arg = (s * (s - a) * (s - b) * (s - c));

% Make sure we don't roundoff to a negative
if(arg < 0) 
   fprintf('STOPPING IN HERON \n');
   [a b c s arg]
   stop
end

area = sqrt(arg);


