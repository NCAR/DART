function quadratic_eq_coeffs

% See also ../doc/{roots_of_m_equation.pptx,cubed_sphere_algorithm.pptx}

logfid = fopen('bb','wt');

% Matrix from http://www.particleincell.com/blog/2012/quad-interpolation/
% describing mapping from an arbitrary quadrilateral to the unit square.

A=[1 0 0 0;1 1 0 0;1 1 1 1;1 0 1 0];
AI = inv(A);

% Explore some of the possible shapes of quads.
for (x2=1:2:10)
for (x3=1:2:10)
for (x4=1:2:10)
for (y2=1:2:10)
for (y3=1:2:10)
for (y4=1:2:10)
   % Make sure that quad is not concave using slopes to the 3 corners.
   if ((y4/x4 >= y3/x3) && (y3/x3 >= y2/x2))
      px = [0, x2, x3, x4];
      py = [0, y2, y3, y4];
      a = AI*px';
      b = AI*py';

      for (o=1:4)
         % Ob locations at the corners of the quad
         x=px(o);
         y=py(o);
         % Put corner 1 at the origin -> a(1) = 0 and b(1) = 0.
         % Simplified bb equation:
         bb =  a(2)*b(3) - a(3)*b(2) + x*b(4) - y*a(4);
         cc =  x*b(2) - y*a(2);

         if (bb < 0 || cc*bb > 0)
            string1 = sprintf('%-g %-g   %-g %-g   %-g %-g  (ox,oy)=(%-g %-g)  bb=%-g   cc=%-g' ...
                              ,x2,y2,x3,y3,x4,y4,x,y,bb,cc);
            fprintf (logfid, '%s \n',string1);
         end
      end
   end
end
end
end
end
end
end
 
fclose(logfid);
=================================================================================

Corner 1 is the origin (0 0)

Note that x2=y2 and x3=y3  means that 2 sides are colinear,
   which is a triangle, rather than a quadrilateral.

Also, x1=y1=1 is (much) closer to corner 1 (the origin)
   than corners 3 and 4 are to corner 2 in many examples.  
   So the quad has a very short side; again distorting towards a triangle.

Corner(x y)      candidate point
2     3     4    (x,y) 
---   ---   ---  -------------  -----   -----
1 1   3 3   1 3  (ox,oy)=(1 3)  bb=-2   cc=-2 
1 1   3 3   1 5  (ox,oy)=(1 5)  bb=-4   cc=-4 
1 1   3 3   1 7  (ox,oy)=(1 7)  bb=-6   cc=-6 
1 1   3 3   1 9  (ox,oy)=(1 9)  bb=-8   cc=-8 
1 1   3 5   1 5  (ox,oy)=(3 5)  bb=-4   cc=-2 
1 1   3 5   1 5  (ox,oy)=(1 5)  bb=-2   cc=-4 
1 1   3 5   1 7  (ox,oy)=(3 5)  bb=-8   cc=-2 
1 1   3 5   1 7  (ox,oy)=(1 7)  bb=-4   cc=-6 
1 1   3 5   1 9  (ox,oy)=(3 5)  bb=-12   cc=-2 
1 1   3 5   1 9  (ox,oy)=(1 9)  bb=-6   cc=-8 
1 1   3 7   1 7  (ox,oy)=(3 7)  bb=-4   cc=-4 
1 1   3 7   1 7  (ox,oy)=(1 7)  bb=-2   cc=-6 
1 1   3 7   1 9  (ox,oy)=(3 7)  bb=-8   cc=-4 
1 1   3 7   1 9  (ox,oy)=(1 9)  bb=-4   cc=-8 
1 1   3 9   1 9  (ox,oy)=(3 9)  bb=-4   cc=-6 
1 1   3 9   1 9  (ox,oy)=(1 9)  bb=-2   cc=-8 
1 3   3 9   1 5  (ox,oy)=(1 5)  bb=-2   cc=-2 
1 3   3 9   1 7  (ox,oy)=(1 7)  bb=-4   cc=-4 
1 3   3 9   1 9  (ox,oy)=(1 9)  bb=-6   cc=-6 
1 1   3 3   3 5  (ox,oy)=(3 5)  bb=-2   cc=-2 
1 1   3 3   3 7  (ox,oy)=(3 7)  bb=-4   cc=-4 
1 1   3 3   3 9  (ox,oy)=(3 9)  bb=-6   cc=-6 
1 1   3 5   3 9  (ox,oy)=(3 5)  bb=-4   cc=-2 
1 1   3 3   5 7  (ox,oy)=(5 7)  bb=-2   cc=-2 
1 1   3 3   5 9  (ox,oy)=(5 9)  bb=-4   cc=-4 
1 1   3 3   7 9  (ox,oy)=(7 9)  bb=-2   cc=-2 
1 1   5 5   1 3  (ox,oy)=(1 3)  bb=-6   cc=-2 
1 1   5 5   1 5  (ox,oy)=(1 5)  bb=-12   cc=-4 
1 1   5 5   1 7  (ox,oy)=(1 7)  bb=-18   cc=-6 
1 1   5 5   1 9  (ox,oy)=(1 9)  bb=-24   cc=-8 
1 1   5 7   1 3  (ox,oy)=(5 7)  bb=-4   cc=-2 
1 1   5 7   1 3  (ox,oy)=(1 3)  bb=-4   cc=-2 
1 1   5 7   1 5  (ox,oy)=(5 7)  bb=-12   cc=-2 
1 1   5 7   1 5  (ox,oy)=(1 5)  bb=-10   cc=-4 
1 1   5 7   1 7  (ox,oy)=(5 7)  bb=-20   cc=-2 
1 1   5 7   1 7  (ox,oy)=(1 7)  bb=-16   cc=-6 
1 1   5 7   1 9  (ox,oy)=(5 7)  bb=-28   cc=-2 
1 1   5 7   1 9  (ox,oy)=(1 9)  bb=-22   cc=-8 
1 1   5 9   1 3  (ox,oy)=(1 3)  bb=-2   cc=-2 
1 1   5 9   1 5  (ox,oy)=(5 9)  bb=-8   cc=-4 
1 1   5 9   1 5  (ox,oy)=(1 5)  bb=-8   cc=-4 
1 1   5 9   1 7  (ox,oy)=(5 9)  bb=-16   cc=-4 
1 1   5 9   1 7  (ox,oy)=(1 7)  bb=-14   cc=-6 
1 1   5 9   1 9  (ox,oy)=(5 9)  bb=-24   cc=-4 
1 1   5 9   1 9  (ox,oy)=(1 9)  bb=-20   cc=-8 
1 1   5 5   3 5  (ox,oy)=(3 5)  bb=-6   cc=-2 
1 1   5 5   3 7  (ox,oy)=(3 7)  bb=-12   cc=-4 
1 1   5 5   3 9  (ox,oy)=(3 9)  bb=-18   cc=-6 
1 1   5 7   3 7  (ox,oy)=(5 7)  bb=-8   cc=-2 
1 1   5 7   3 7  (ox,oy)=(3 7)  bb=-6   cc=-4 
1 1   5 7   3 9  (ox,oy)=(5 7)  bb=-16   cc=-2 
1 1   5 7   3 9  (ox,oy)=(3 9)  bb=-12   cc=-6 
1 1   5 9   3 9  (ox,oy)=(5 9)  bb=-8   cc=-4 
1 1   5 9   3 9  (ox,oy)=(3 9)  bb=-6   cc=-6 
1 1   5 5   5 7  (ox,oy)=(5 7)  bb=-6   cc=-2 
1 1   5 5   5 9  (ox,oy)=(5 9)  bb=-12   cc=-4 
1 1   5 7   5 9  (ox,oy)=(5 7)  bb=-4   cc=-2 
1 1   5 7   5 9  (ox,oy)=(5 9)  bb=-2   cc=-4 
1 1   5 5   7 9  (ox,oy)=(7 9)  bb=-6   cc=-2 
1 1   7 7   1 3  (ox,oy)=(1 3)  bb=-10   cc=-2 
1 1   7 7   1 5  (ox,oy)=(1 5)  bb=-20   cc=-4 
1 1   7 7   1 7  (ox,oy)=(1 7)  bb=-30   cc=-6 
1 1   7 7   1 9  (ox,oy)=(1 9)  bb=-40   cc=-8 
1 1   7 9   1 3  (ox,oy)=(7 9)  bb=-8   cc=-2 
1 1   7 9   1 3  (ox,oy)=(1 3)  bb=-8   cc=-2 
1 1   7 9   1 5  (ox,oy)=(7 9)  bb=-20   cc=-2 
1 1   7 9   1 5  (ox,oy)=(1 5)  bb=-18   cc=-4 
1 1   7 9   1 7  (ox,oy)=(7 9)  bb=-32   cc=-2 
1 1   7 9   1 7  (ox,oy)=(1 7)  bb=-28   cc=-6 
1 1   7 9   1 9  (ox,oy)=(7 9)  bb=-44   cc=-2 
1 1   7 9   1 9  (ox,oy)=(1 9)  bb=-38   cc=-8 
1 1   7 7   3 5  (ox,oy)=(3 5)  bb=-10   cc=-2 
1 1   7 7   3 7  (ox,oy)=(3 7)  bb=-20   cc=-4 
1 1   7 7   3 9  (ox,oy)=(3 9)  bb=-30   cc=-6 
1 1   7 9   3 5  (ox,oy)=(7 9)  bb=-4   cc=-2 
1 1   7 9   3 5  (ox,oy)=(3 5)  bb=-4   cc=-2 
1 1   7 9   3 7  (ox,oy)=(7 9)  bb=-16   cc=-2 
1 1   7 9   3 7  (ox,oy)=(3 7)  bb=-14   cc=-4 
1 1   7 9   3 9  (ox,oy)=(7 9)  bb=-28   cc=-2 
1 1   7 9   3 9  (ox,oy)=(3 9)  bb=-24   cc=-6 
1 1   7 7   5 7  (ox,oy)=(5 7)  bb=-10   cc=-2 
1 1   7 7   5 9  (ox,oy)=(5 9)  bb=-20   cc=-4 
1 1   7 9   5 9  (ox,oy)=(7 9)  bb=-12   cc=-2 
1 1   7 9   5 9  (ox,oy)=(5 9)  bb=-10   cc=-4 
1 1   7 7   7 9  (ox,oy)=(7 9)  bb=-10   cc=-2 
1 1   9 9   1 3  (ox,oy)=(1 3)  bb=-14   cc=-2 
1 1   9 9   1 5  (ox,oy)=(1 5)  bb=-28   cc=-4 
1 1   9 9   1 7  (ox,oy)=(1 7)  bb=-42   cc=-6 
1 1   9 9   1 9  (ox,oy)=(1 9)  bb=-56   cc=-8 
1 1   9 9   3 5  (ox,oy)=(3 5)  bb=-14   cc=-2 
1 1   9 9   3 7  (ox,oy)=(3 7)  bb=-28   cc=-4 
1 1   9 9   3 9  (ox,oy)=(3 9)  bb=-42   cc=-6 
1 1   9 9   5 7  (ox,oy)=(5 7)  bb=-14   cc=-2 
1 1   9 9   5 9  (ox,oy)=(5 9)  bb=-28   cc=-4 
1 1   9 9   7 9  (ox,oy)=(7 9)  bb=-14   cc=-2 
3 1   5 3   1 5  (ox,oy)=(5 3)  bb=-4   cc=-4 
3 1   5 3   1 7  (ox,oy)=(5 3)  bb=-8   cc=-4 
3 1   5 3   1 9  (ox,oy)=(5 3)  bb=-12   cc=-4 
3 1   5 5   1 9  (ox,oy)=(5 5)  bb=-4   cc=-10 
3 3   5 7   1 7  (ox,oy)=(5 7)  bb=-4   cc=-6 
3 3   5 7   1 9  (ox,oy)=(5 7)  bb=-8   cc=-6 
3 5   5 9   1 5  (ox,oy)=(5 9)  bb=-4   cc=-2 
3 5   5 9   1 7  (ox,oy)=(5 9)  bb=-8   cc=-2 
3 5   5 9   1 9  (ox,oy)=(5 9)  bb=-12   cc=-2 
3 1   5 3   3 7  (ox,oy)=(5 3)  bb=-4   cc=-4 
3 1   5 3   3 9  (ox,oy)=(5 3)  bb=-8   cc=-4 
3 5   5 9   3 9  (ox,oy)=(5 9)  bb=-4   cc=-2 
3 1   5 3   5 9  (ox,oy)=(5 3)  bb=-4   cc=-4 
3 1   7 3   1 3  (ox,oy)=(7 3)  bb=-8   cc=-2 
3 1   7 3   1 3  (ox,oy)=(1 3)  bb=-2   cc=-8 
3 1   7 3   1 5  (ox,oy)=(7 3)  bb=-16   cc=-2 
3 1   7 3   1 5  (ox,oy)=(1 5)  bb=-4   cc=-14 
3 1   7 3   1 7  (ox,oy)=(7 3)  bb=-24   cc=-2 
3 1   7 3   1 7  (ox,oy)=(1 7)  bb=-6   cc=-20 
3 1   7 3   1 9  (ox,oy)=(7 3)  bb=-32   cc=-2 
3 1   7 3   1 9  (ox,oy)=(1 9)  bb=-8   cc=-26 
3 1   7 5   1 5  (ox,oy)=(7 5)  bb=-8   cc=-8 
3 1   7 5   1 5  (ox,oy)=(1 5)  bb=-2   cc=-14 
3 1   7 5   1 7  (ox,oy)=(7 5)  bb=-16   cc=-8 
3 1   7 5   1 7  (ox,oy)=(1 7)  bb=-4   cc=-20 
3 1   7 5   1 9  (ox,oy)=(7 5)  bb=-24   cc=-8 
3 1   7 5   1 9  (ox,oy)=(1 9)  bb=-6   cc=-26 
3 1   7 7   1 7  (ox,oy)=(7 7)  bb=-8   cc=-14 
3 1   7 7   1 7  (ox,oy)=(1 7)  bb=-2   cc=-20 
3 1   7 7   1 9  (ox,oy)=(7 7)  bb=-16   cc=-14 
3 1   7 7   1 9  (ox,oy)=(1 9)  bb=-4   cc=-26 
3 1   7 9   1 9  (ox,oy)=(7 9)  bb=-8   cc=-20 
3 1   7 9   1 9  (ox,oy)=(1 9)  bb=-2   cc=-26 
3 3   7 7   1 3  (ox,oy)=(1 3)  bb=-2   cc=-6 
3 3   7 7   1 5  (ox,oy)=(1 5)  bb=-4   cc=-12 
3 3   7 7   1 7  (ox,oy)=(1 7)  bb=-6   cc=-18 
3 3   7 7   1 9  (ox,oy)=(1 9)  bb=-8   cc=-24 
3 3   7 9   1 5  (ox,oy)=(7 9)  bb=-8   cc=-6 
3 3   7 9   1 5  (ox,oy)=(1 5)  bb=-2   cc=-12 
3 3   7 9   1 7  (ox,oy)=(7 9)  bb=-16   cc=-6 
3 3   7 9   1 7  (ox,oy)=(1 7)  bb=-4   cc=-18 
3 3   7 9   1 9  (ox,oy)=(7 9)  bb=-24   cc=-6 
3 3   7 9   1 9  (ox,oy)=(1 9)  bb=-6   cc=-24 
3 1   7 3   3 3  (ox,oy)=(7 3)  bb=-4   cc=-2 
3 1   7 3   3 5  (ox,oy)=(7 3)  bb=-12   cc=-2 
3 1   7 3   3 5  (ox,oy)=(3 5)  bb=-2   cc=-12 
3 1   7 3   3 7  (ox,oy)=(7 3)  bb=-20   cc=-2 
3 1   7 3   3 7  (ox,oy)=(3 7)  bb=-4   cc=-18 
3 1   7 3   3 9  (ox,oy)=(7 3)  bb=-28   cc=-2 
3 1   7 3   3 9  (ox,oy)=(3 9)  bb=-6   cc=-24 
3 1   7 5   3 7  (ox,oy)=(7 5)  bb=-8   cc=-8 
3 1   7 5   3 9  (ox,oy)=(7 5)  bb=-16   cc=-8 
3 1   7 7   3 9  (ox,oy)=(7 7)  bb=-4   cc=-14 
3 3   7 7   3 5  (ox,oy)=(3 5)  bb=-2   cc=-6 
3 3   7 7   3 7  (ox,oy)=(3 7)  bb=-4   cc=-12 
3 3   7 7   3 9  (ox,oy)=(3 9)  bb=-6   cc=-18 
3 3   7 9   3 7  (ox,oy)=(7 9)  bb=-4   cc=-6 
3 3   7 9   3 9  (ox,oy)=(7 9)  bb=-12   cc=-6 
3 1   7 3   5 5  (ox,oy)=(7 3)  bb=-8   cc=-2 
3 1   7 3   5 7  (ox,oy)=(7 3)  bb=-16   cc=-2 
3 1   7 3   5 7  (ox,oy)=(5 7)  bb=-2   cc=-16 
3 1   7 3   5 9  (ox,oy)=(7 3)  bb=-24   cc=-2 
3 1   7 3   5 9  (ox,oy)=(5 9)  bb=-4   cc=-22 
3 1   7 5   5 9  (ox,oy)=(7 5)  bb=-8   cc=-8 
3 3   7 7   5 7  (ox,oy)=(5 7)  bb=-2   cc=-6 
3 3   7 7   5 9  (ox,oy)=(5 9)  bb=-4   cc=-12 
3 1   7 3   7 5  (ox,oy)=(7 3)  bb=-4   cc=-2 
3 1   7 3   7 7  (ox,oy)=(7 3)  bb=-12   cc=-2 
3 1   7 3   7 9  (ox,oy)=(7 3)  bb=-20   cc=-2 
3 1   7 3   7 9  (ox,oy)=(7 9)  bb=-2   cc=-20 
3 3   7 7   7 9  (ox,oy)=(7 9)  bb=-2   cc=-6 
3 1   7 3   9 7  (ox,oy)=(7 3)  bb=-8   cc=-2 
3 1   7 3   9 9  (ox,oy)=(7 3)  bb=-16   cc=-2 
3 1   9 3   1 1  (ox,oy)=(1 1)  bb=-2   cc=-2 
3 1   9 3   1 3  (ox,oy)=(1 3)  bb=-8   cc=-8 
3 1   9 3   1 5  (ox,oy)=(1 5)  bb=-14   cc=-14 
3 1   9 3   1 7  (ox,oy)=(1 7)  bb=-20   cc=-20 
3 1   9 3   1 9  (ox,oy)=(1 9)  bb=-26   cc=-26 
3 1   9 5   1 3  (ox,oy)=(9 5)  bb=-8   cc=-6 
3 1   9 5   1 3  (ox,oy)=(1 3)  bb=-6   cc=-8 
3 1   9 5   1 5  (ox,oy)=(9 5)  bb=-20   cc=-6 
3 1   9 5   1 5  (ox,oy)=(1 5)  bb=-12   cc=-14 
3 1   9 5   1 7  (ox,oy)=(9 5)  bb=-32   cc=-6 
3 1   9 5   1 7  (ox,oy)=(1 7)  bb=-18   cc=-20 
3 1   9 5   1 9  (ox,oy)=(9 5)  bb=-44   cc=-6 
3 1   9 5   1 9  (ox,oy)=(1 9)  bb=-24   cc=-26 
3 1   9 7   1 3  (ox,oy)=(1 3)  bb=-4   cc=-8 
3 1   9 7   1 5  (ox,oy)=(9 7)  bb=-12   cc=-12 
3 1   9 7   1 5  (ox,oy)=(1 5)  bb=-10   cc=-14 
3 1   9 7   1 7  (ox,oy)=(9 7)  bb=-24   cc=-12 
3 1   9 7   1 7  (ox,oy)=(1 7)  bb=-16   cc=-20 
3 1   9 7   1 9  (ox,oy)=(9 7)  bb=-36   cc=-12 
3 1   9 7   1 9  (ox,oy)=(1 9)  bb=-22   cc=-26 
3 1   9 9   1 3  (ox,oy)=(1 3)  bb=-2   cc=-8 
3 1   9 9   1 5  (ox,oy)=(9 9)  bb=-4   cc=-18 
3 1   9 9   1 5  (ox,oy)=(1 5)  bb=-8   cc=-14 
3 1   9 9   1 7  (ox,oy)=(9 9)  bb=-16   cc=-18 
3 1   9 9   1 7  (ox,oy)=(1 7)  bb=-14   cc=-20 
3 1   9 9   1 9  (ox,oy)=(9 9)  bb=-28   cc=-18 
3 1   9 9   1 9  (ox,oy)=(1 9)  bb=-20   cc=-26 
3 3   9 9   1 3  (ox,oy)=(1 3)  bb=-6   cc=-6 
3 3   9 9   1 5  (ox,oy)=(1 5)  bb=-12   cc=-12 
3 3   9 9   1 7  (ox,oy)=(1 7)  bb=-18   cc=-18 
3 3   9 9   1 9  (ox,oy)=(1 9)  bb=-24   cc=-24 
3 1   9 3   3 3  (ox,oy)=(3 3)  bb=-6   cc=-6 
3 1   9 3   3 5  (ox,oy)=(3 5)  bb=-12   cc=-12 
3 1   9 3   3 7  (ox,oy)=(3 7)  bb=-18   cc=-18 
3 1   9 3   3 9  (ox,oy)=(3 9)  bb=-24   cc=-24 
3 1   9 5   3 5  (ox,oy)=(9 5)  bb=-12   cc=-6 
3 1   9 5   3 5  (ox,oy)=(3 5)  bb=-6   cc=-12 
3 1   9 5   3 7  (ox,oy)=(9 5)  bb=-24   cc=-6 
3 1   9 5   3 7  (ox,oy)=(3 7)  bb=-12   cc=-18 
3 1   9 5   3 9  (ox,oy)=(9 5)  bb=-36   cc=-6 
3 1   9 5   3 9  (ox,oy)=(3 9)  bb=-18   cc=-24 
3 1   9 7   3 7  (ox,oy)=(9 7)  bb=-12   cc=-12 
3 1   9 7   3 7  (ox,oy)=(3 7)  bb=-6   cc=-18 
3 1   9 7   3 9  (ox,oy)=(9 7)  bb=-24   cc=-12 
3 1   9 7   3 9  (ox,oy)=(3 9)  bb=-12   cc=-24 
3 1   9 9   3 9  (ox,oy)=(9 9)  bb=-12   cc=-18 
3 1   9 9   3 9  (ox,oy)=(3 9)  bb=-6   cc=-24 
3 3   9 9   3 5  (ox,oy)=(3 5)  bb=-6   cc=-6 
3 3   9 9   3 7  (ox,oy)=(3 7)  bb=-12   cc=-12 
3 3   9 9   3 9  (ox,oy)=(3 9)  bb=-18   cc=-18 
3 1   9 3   5 3  (ox,oy)=(5 3)  bb=-4   cc=-4 
3 1   9 3   5 5  (ox,oy)=(5 5)  bb=-10   cc=-10 
3 1   9 3   5 7  (ox,oy)=(5 7)  bb=-16   cc=-16 
3 1   9 3   5 9  (ox,oy)=(5 9)  bb=-22   cc=-22 
3 1   9 5   5 5  (ox,oy)=(9 5)  bb=-4   cc=-6 
3 1   9 5   5 7  (ox,oy)=(9 5)  bb=-16   cc=-6 
3 1   9 5   5 7  (ox,oy)=(5 7)  bb=-6   cc=-16 
3 1   9 5   5 9  (ox,oy)=(9 5)  bb=-28   cc=-6 
3 1   9 5   5 9  (ox,oy)=(5 9)  bb=-12   cc=-22 
3 1   9 7   5 9  (ox,oy)=(9 7)  bb=-12   cc=-12 
3 1   9 7   5 9  (ox,oy)=(5 9)  bb=-2   cc=-22 
3 3   9 9   5 7  (ox,oy)=(5 7)  bb=-6   cc=-6 
3 3   9 9   5 9  (ox,oy)=(5 9)  bb=-12   cc=-12 
3 1   9 3   7 3  (ox,oy)=(7 3)  bb=-2   cc=-2 
3 1   9 3   7 5  (ox,oy)=(7 5)  bb=-8   cc=-8 
3 1   9 3   7 7  (ox,oy)=(7 7)  bb=-14   cc=-14 
3 1   9 3   7 9  (ox,oy)=(7 9)  bb=-20   cc=-20 
3 1   9 5   7 7  (ox,oy)=(9 5)  bb=-8   cc=-6 
3 1   9 5   7 9  (ox,oy)=(9 5)  bb=-20   cc=-6 
3 1   9 5   7 9  (ox,oy)=(7 9)  bb=-6   cc=-20 
3 3   9 9   7 9  (ox,oy)=(7 9)  bb=-6   cc=-6 
3 1   9 3   9 5  (ox,oy)=(9 5)  bb=-6   cc=-6 
3 1   9 3   9 7  (ox,oy)=(9 7)  bb=-12   cc=-12 
3 1   9 3   9 9  (ox,oy)=(9 9)  bb=-18   cc=-18 
3 1   9 5   9 9  (ox,oy)=(9 5)  bb=-12   cc=-6 
5 1   7 3   1 7  (ox,oy)=(7 3)  bb=-4   cc=-8 
5 1   7 3   1 9  (ox,oy)=(7 3)  bb=-8   cc=-8 
5 3   7 5   1 5  (ox,oy)=(7 5)  bb=-4   cc=-4 
5 3   7 5   1 7  (ox,oy)=(7 5)  bb=-8   cc=-4 
5 3   7 5   1 9  (ox,oy)=(7 5)  bb=-12   cc=-4 
5 5   7 9   1 9  (ox,oy)=(7 9)  bb=-4   cc=-10 
5 1   7 3   3 9  (ox,oy)=(7 3)  bb=-4   cc=-8 
5 3   7 5   3 7  (ox,oy)=(7 5)  bb=-4   cc=-4 
5 3   7 5   3 9  (ox,oy)=(7 5)  bb=-8   cc=-4 
5 3   7 5   5 9  (ox,oy)=(7 5)  bb=-4   cc=-4 
5 1   9 3   1 3  (ox,oy)=(9 3)  bb=-4   cc=-6 
5 1   9 3   1 5  (ox,oy)=(9 3)  bb=-12   cc=-6 
5 1   9 3   1 7  (ox,oy)=(9 3)  bb=-20   cc=-6 
5 1   9 3   1 9  (ox,oy)=(9 3)  bb=-28   cc=-6 
5 1   9 5   1 7  (ox,oy)=(9 5)  bb=-8   cc=-16 
5 1   9 5   1 9  (ox,oy)=(9 5)  bb=-16   cc=-16 
5 1   9 7   1 9  (ox,oy)=(9 7)  bb=-4   cc=-26 
5 3   9 7   1 5  (ox,oy)=(9 7)  bb=-8   cc=-8 
5 3   9 7   1 7  (ox,oy)=(9 7)  bb=-16   cc=-8 
5 3   9 7   1 9  (ox,oy)=(9 7)  bb=-24   cc=-8 
5 3   9 9   1 7  (ox,oy)=(9 9)  bb=-4   cc=-18 
5 3   9 9   1 9  (ox,oy)=(9 9)  bb=-12   cc=-18 
5 1   9 3   3 5  (ox,oy)=(9 3)  bb=-8   cc=-6 
5 1   9 3   3 7  (ox,oy)=(9 3)  bb=-16   cc=-6 
5 1   9 3   3 9  (ox,oy)=(9 3)  bb=-24   cc=-6 
5 1   9 5   3 9  (ox,oy)=(9 5)  bb=-8   cc=-16 
5 3   9 7   3 7  (ox,oy)=(9 7)  bb=-8   cc=-8 
5 3   9 7   3 9  (ox,oy)=(9 7)  bb=-16   cc=-8 
5 1   9 3   5 5  (ox,oy)=(9 3)  bb=-4   cc=-6 
5 1   9 3   5 7  (ox,oy)=(9 3)  bb=-12   cc=-6 
5 1   9 3   5 9  (ox,oy)=(9 3)  bb=-20   cc=-6 
5 3   9 7   5 9  (ox,oy)=(9 7)  bb=-8   cc=-8 
5 1   9 3   7 7  (ox,oy)=(9 3)  bb=-8   cc=-6 
5 1   9 3   7 9  (ox,oy)=(9 3)  bb=-16   cc=-6 
5 1   9 3   9 7  (ox,oy)=(9 3)  bb=-4   cc=-6 
5 1   9 3   9 9  (ox,oy)=(9 3)  bb=-12   cc=-6 
7 1   9 3   1 9  (ox,oy)=(9 3)  bb=-4   cc=-12 
7 3   9 5   1 7  (ox,oy)=(9 5)  bb=-4   cc=-8 
7 3   9 5   1 9  (ox,oy)=(9 5)  bb=-8   cc=-8 
7 5   9 7   1 5  (ox,oy)=(9 7)  bb=-4   cc=-4 
7 5   9 7   1 7  (ox,oy)=(9 7)  bb=-8   cc=-4 
7 5   9 7   1 9  (ox,oy)=(9 7)  bb=-12   cc=-4 
7 3   9 5   3 9  (ox,oy)=(9 5)  bb=-4   cc=-8 
7 5   9 7   3 7  (ox,oy)=(9 7)  bb=-4   cc=-4 
7 5   9 7   3 9  (ox,oy)=(9 7)  bb=-8   cc=-4 
7 5   9 7   5 9  (ox,oy)=(9 7)  bb=-4   cc=-4 
