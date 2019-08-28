function x = norm_inv(p)
%   Computes inverse CDF for the normal distribution,
%   evaluated at either the scalar value or array P.
%
% based on the code in DART assim_tools_mod.f90, which is
% in turn based on http://home.online.no/~pjacklam/notes/invnorm

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% constants where precision is vital
a1 = -39.69683028665376;
a2 =  220.9460984245205;
a3 = -275.9285104469687;
a4 =  138.357751867269;
a5 = -30.66479806614716;
a6 =  2.506628277459239;
b1 = -54.4760987982241;
b2 =  161.5858368580409;
b3 = -155.6989798598866;
b4 =  66.80131188771972;
b5 = -13.28068155288572;
c1 = -0.007784894002430293;
c2 = -0.3223964580411365;
c3 = -2.400758277161838;
c4 = -2.549732539343734;
c5 =  4.374664141464968;
c6 =  2.938163982698783;
d1 =  0.007784695709041462;
d2 =  0.3224671290700398;
d3 =  2.445134137142996;
d4 =  3.754408661907416;

% Split into an inner and two outer regions which have separate fits
p_low  = 0.02425;
p_high = 1 - p_low;

% if P is scalar, simplify
if (numel(p) == 1)
 if (p < p_low)
  q = sqrt(-2.0 * log(p));
  x = (((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) / ...
       ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0);
 elseif (p >= p_low) && (p <= p_high)
  q = p - 0.5;
  r = q*q;
  x = (((((a1*r + a2)*r + a3)*r + a4)*r + a5)*r + a6)*q / ...
      (((((b1*r + b2)*r + b3)*r + b4)*r + b5)*r + 1.0);

 else % (p > p_high)
  q = sqrt(-2.0 * log(1.0 - p));
  x = -(((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) / ...
        ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0);
 end

else % if P is an array, so do it element by element

 inds = (p < p_low);
 q(inds) = sqrt(-2.0 * log(p(inds)));
 x(inds) = (((((c1*q(inds) + c2).*q(inds) + c3).*q(inds) + c4).*q(inds) + c5).*q(inds) + c6) ./ ...
            ((((d1*q(inds) + d2).*q(inds) + d3).*q(inds) + d4).*q(inds) + 1.0);

 inds = (p >= p_low) & (p <= p_high);
 q(inds) = p(inds) - 0.5;
 r(inds) = q(inds).*q(inds);
 x(inds) = (((((a1*r(inds) + a2).*r(inds) + a3).*r(inds) + a4).*r(inds) + a5).*r(inds) + a6).*q(inds) ./ ...
           (((((b1*r(inds) + b2).*r(inds) + b3).*r(inds) + b4).*r(inds) + b5).*r(inds) + 1.0);

 inds = (p > p_high);
 q(inds) = sqrt(-2.0 * log(1.0 - p(inds)));
 x(inds) = -(((((c1*q(inds) + c2).*q(inds) + c3).*q(inds) + c4).*q(inds) + c5).*q(inds) + c6) ./ ...
             ((((d1*q(inds) + d2).*q(inds) + d3).*q(inds) + d4).*q(inds) + 1.0);
end

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
