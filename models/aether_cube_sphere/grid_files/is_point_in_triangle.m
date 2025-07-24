% Determines if the projection of a point p onto the plane of a triangle with vertices
% v1, v2 and v3 is inside the triangle or not. Computes the areas of each of the triangles
% between p and a pair of vertices. These should sum to the area of the triangle if p
% is inside and be larger than that if p is outside. 

function [inside] = is_point_in_triangle(v1, v2, v3, p)


% Get the projection of the point p onto the plane containing the triangle
% Start by getting perpendicular vector to plane by cross product
perp = cross((v1-v2), v2-v3);
% Get unit vector in direction of perp
unit_perp = perp / norm(perp);
% Projection of vector from v1 to p on the unit perp vector is how much to move to get to plane
offset = dot((p-v1), unit_perp);
p_proj = p - offset * unit_perp;

% Compute lengths of the sides
len_s1 = norm((v1 - v2));
len_s2 = norm((v3 - v2));
len_s3 = norm((v1 - v3));

% Compute the lengths from the point p
len_p1 = norm((p_proj - v1));
len_p2 = norm((p_proj - v2));
len_p3 = norm((p_proj - v3));

% Area of triangle
at = heron(len_s1, len_s2, len_s3);

% Compute areas of sub triangles
at1 = heron(len_p1, len_p2, len_s1);
at2 = heron(len_p2, len_p3, len_s2);
at3 = heron(len_p3, len_p1, len_s3);

% Difference between sub triangles and the triangle area
area_dif = at1 + at2 + at3 - at;

% Quadrilaterals on the interior of the cube sphere sides are really spherical quads,
% There sides are great circles. This routine assumes that the triangles composing the quads
% have straight sides in regular space. The algorithm finds points that are inside the 
% spherical quads. These quads actually 'bulge' out compared to the regular sides, so it is possible
% to have points that are inside the spherical quad but just barely outside of the regular
% quads. This threshold is tuned so that these points still show as inside. The tuning is for
% np = 18 (number of points along a grid face is 18). Fewer points might require a larger
% threshold while more points might be okay with a smaller one.
threshold = 0.002;

dif_frac = area_dif / at;
inside = abs(dif_frac) < threshold;

return

