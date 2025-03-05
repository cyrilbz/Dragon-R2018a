function s = multiOrientSkel(BW, msk, prc_thresh)
% MULTIORIENTSKEL  Generates skeleton in masked region from greyscale image
%                  Image is rotated so that the skeleton lies along x-axis
%
%   @input: BW         - greyscale image to be skeletonised
%           msk        - binary mask matching dimensionality of im
%           prc_thresh - percentile threshold passed to getBinary
%
%   @output: s - skeletonised image of im, rotated to lie along x-axis
%
%   Skeletonisation is performed with and without rotation and the results
%   combined, mitigating data loss from rotating discrete data

    p = getParams();

    ang = getOrientation(msk);
    ang_rt = rad2deg(-ang);
    
    im_rt = imrotate(BW,ang_rt);
    msk_rt = imrotate(msk,ang_rt);
    
    bin = getBinary(BW, msk, prc_thresh);
    bin_rt = getBinary(im_rt, msk_rt, prc_thresh);

    skel_rt = bwskel(bin_rt, 'MinBranchLength', p.minBranchLength);
    skel = bwskel(bin, 'MinBranchLength', p.minBranchLength);
    
    skel2 = imrotate(skel,ang_rt);
    s = bwskel((skel2 | skel_rt) > 0, 'MinBranchLength', p.minBranchLength);
end