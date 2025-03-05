function bp = getBranchpoints(im_skel)
% GETBRANCHPOINTS  Determines branchpoints (with refinement over BWMORPH)
%
%   @input: im_skel - a binary skeleton
%
%   @output: bp - location of all the branchpoints (binary image)

    im_branch = bwmorph(im_skel, 'branchpoints');
    im_skel_seg = im_skel - im_branch;      %remove branchpoints from skeleton
    im_branch2 = bwmorph(im_skel_seg, 'branchpoints');
    im_branch_sum = im_branch + im_branch2;
    pairs = bwmorph(im_branch_sum,'clean');
    bp = xor(im_branch,pairs);
end