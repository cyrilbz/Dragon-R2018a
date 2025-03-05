function curve = measureCurve(a,b,c)
% MEASURECURVE Determines the absolute value of the curvature of three
% points. The Menger curvature is defined here: 
% https://en.wikipedia.org/wiki/Menger_curvature#Definition
%
%   @input: a,b,c - [x,y] coordinates of the 3 points
%
%   @output: curve - Curvature from those three points

    Area = triArea(a, b, c);
    if(Area == 0)  % colinear points, R->inf, c->0
        curve = 0;
    else
        %curve = (4*Area)/abs(pdist2(a,b)*pdist2(b,c)*pdist2(c,a));
        curve = (4*Area)/abs(norm(a-b)*norm(b-c)*norm(c-a));
    end
end