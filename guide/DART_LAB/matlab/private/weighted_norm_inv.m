
function [x] = weighted_norm_inv(alpha, mean, sd, p)

% Find the value of x for which the cdf of a N(mean, sd)
% multiplied times alpha has value p.

% Can search in a standard normal, then multiply by sd at end and add mean
% Divide p by alpha to get the right place for weighted normal
np = p / alpha;

% Find spot in standard normal
x = norm_inv(np);

% Add in the mean and normalize by sd
x = mean + x * sd;


