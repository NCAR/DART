function [y] = donorm(x, mu, sigma)
%  computes a gaussian (normal) probabilty distribution function
%  for the points of X with a given mean (mu) and standard deviation (sigma)
% 
% normal plot, y given x:
%  y = (1 / sigma * sqrt(2*pi)) * e ^^ ((-1/2 * ((x-mu) / sigma)^^2)
% or
%  g(x) = \frac{1}{\sigma\sqrt{2\pi}} e^{ -\frac{1}{2}\left(\frac{x-\mu}{\sigma}\right)^2 }. 

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$


e = 2.71828182845904523536;

basen = (1.0 / (sigma * sqrt(2*pi)));
expon = -0.50 .* (((x-mu) / sigma).^2 );

y = basen .* (e .^ expon);

end

 
