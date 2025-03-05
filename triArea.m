function A = triArea(a,b,c)
% TRIAREA  returns the area between three points that are not colinear
%
%   @input: a,b,c - three coordinates in the format (x,y)
%
%   @output: A - area contained within three points forming a triangle

    A = a(1) * (b(2) - c(2));
    A = A + b(1) * (c(2) - a(2));
    A = A + c(1) * (a(2) - b(2));
    A = abs(A)/2.0;
end