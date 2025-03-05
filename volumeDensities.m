function np = volumeDensities(np, im_stack)
% VOLUMEDENSITIES Scale ROI area by number of z-stack images and step size.
% Calculate VolDensities based on this 3D volume estimate
%
%   @input: np       - current net_props structure (state: end of angleZ)
%           im_stack - full z-stack of images
%
%   @output: np - updated net_props structure
%
%   NOTE-1: Will not update the np object if im_stack is not 3 dimensional
%
%   NOTE-2: The volume is purely extrapolated using the number of z-slices
%   and the distance between them relative to the pixel size. This DOES NOT
%   account for cell topology and gives an approimate volume in pixels
%   cubed.



    if(ndims(im_stack) ~= 3)
        return
    end
    np.z_stack_size = size(im_stack,1);
    np.cellVolume = np.cellSize * (np.z_scale * np.z_stack_size);

    size_to_vol_ratio = np.cellSize / np.cellVolume;
    np.actinVolDen = np.actinDen * size_to_vol_ratio;
    np.branchVolDen = np.branchDen  * size_to_vol_ratio;
    np.skelVolDensity = np.skelDensity * size_to_vol_ratio;

end