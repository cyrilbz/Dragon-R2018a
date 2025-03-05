function stk = grabImageStack(filename)
% GRABIMAGESTACK  This function generates a 3D array from a stack of 2D images
%
%   @input: filename - the name(+path) of the stack of images
%
%   @output: stk - 3D array of 2D image stack. First index goes through stack
%
%   This was tested using .tif files with all slices in a single TIFF image

    info = imfinfo(filename);
    num_imgs = numel(info);
    im_test = imread(filename);
    [imY, imX] = size(im_test);
    type = class(im_test);
    stk = zeros(num_imgs, imY, imX, type);
    for im_no = 1:num_imgs
        stk(im_no,:,:) = imread(filename, im_no, 'Info', info);
    end
end