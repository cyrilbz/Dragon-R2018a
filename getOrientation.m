function ori = getOrientation(I)
% GETORIENTATION  calculates the mean orientation (radians) of the regions
% from the input image
%
%   @input: I - image used to determine region orientation
%
%   @output: ori - Mean orientation of regions (in radians) in I
%
%   Using a binary mask is an effective way to determine the orientation of a cell
%   or a specific region

    imRawStats = regionprops(I, 'orientation');
    len = length(imRawStats);
    sum = 0.0;
    for i = 1:len
        sum = sum + imRawStats(i).Orientation;
    end
    ori = sum/len;
    ori = deg2rad(ori);
end