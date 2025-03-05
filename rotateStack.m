function rt_stk = rotateStack(stk, ang)
% ROTATESTACK  This function rotates every image in an image stack
%
%   @input: stk - the 3D stack of a 2D image series
%           ang - the angle for rotation using imrotate()
%
%   @output: rt_stk - the image stack with all images rotated and resized
%
%   Uses the same angle definition for rotation as IMROTATE

    [stackLen, ~, ~] = size(stk);
    im_test = imrotate(squeeze(stk(1,:,:)), ang);
    [sizeY, sizeX] = size(im_test);
    rt_stk = zeros(stackLen, sizeY, sizeX);
    for i = 1:stackLen
        pre_rt = squeeze(stk(i,:,:));
        post_rt = imrotate(pre_rt, ang);
        rt_stk(i,:,:) = post_rt;
    end
end