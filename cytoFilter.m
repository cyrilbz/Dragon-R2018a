function im_filter = cytoFilter(path)
% CYTOFILTER  Applies a 'rolling ball'-like background subtraction to the image
%
%   @input: path - the path containing the PEN3 and actin image files
%
%   @output: filtered actin image
%
%   mask is saved to the 'path' specified as input

    base_filename = "actin_original.tif";
    filter_filename = "actin_filter.tif";
    
    p = getParams();

    im_raw = imread(char(path + base_filename));
    im1 = to8Bit(im_raw);
    im_filter = imtophat(im1,strel('sphere',p.rb_radius));
    imwrite(im_filter,char( path + filter_filename));
end