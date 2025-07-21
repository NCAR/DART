function [area] = heron(a, b, c)

s = (a + b + c) /2;

arg = (s * (s - a) * (s - b) * (s -c));

% Make sure we don't roundoff to a negative
if(arg < 0) arg = 0; end

area = sqrt(arg);

if(imag(area) ~= 0)
   format long
   [a b c s area]
   stop
end

