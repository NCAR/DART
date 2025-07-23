% Computes Herons formula to get area of triangle from lenghts of sides
% Super accuracy is not needed in the area calculation here

function [area] = heron(a, b, c)

s = (a + b + c) /2;

arg = (s * (s - a) * (s - b) * (s - c));

if(arg < -1e-10)
   fprintf('STOPPING IN HERON FOR NEGATIVE AREA \n');
   [a b c s arg]
   stop
end
   
% Make sure we don't roundoff to a negative
if(arg <= 0) 
   area = 0;
   return
end

area = sqrt(arg);
