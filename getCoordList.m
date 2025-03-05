function lst = getCoordList(BW)
% GETCOORDLIST  Returns list of coordinates of true pixels in binary image
%
%   @input: BW - binary image
%
%   @output: lst - Nx2 array of N coordinates in the format (x,y)

    [y,x] = find(BW);
    lst = [x y];
end