function[inside] = is_point_in_quad(v, p)

% See if the point is inside this quad; it's inside if it's in one or more contained triangles
inside_t(1) = is_point_in_triangle(v(1, :), v(2, :), v(3, :), p);
inside_t(2) = is_point_in_triangle(v(1, :), v(2, :), v(4, :), p);
inside_t(3) = is_point_in_triangle(v(1, :), v(3, :), v(4, :), p);
inside_t(4) = is_point_in_triangle(v(2, :), v(3, :), v(4, :), p);

inside = any(inside_t);
   

