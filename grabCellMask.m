function msk = grabCellMask(msk_filename, mskBG)
% GRABCELLMASK  Loads a previously saved mask if one is available, or
% allows one to be drawn over the image passed to the function
%
%   @input: msk_filename - filename (and path) for the cell mask. If one is
%           found here, it will be loaded. Otherwise, the one made in the 
%           function will be saved to this
%           mskBG - the image used to draw over in mask creation
%
%   @output: msk - binary cell mask


    if(~isfile(msk_filename))
        imshow(mskBG, []);
        msk = createMask(drawfreehand);
        msk = msk > 0;
        imwrite(msk,path + mask_filename);
        close;
    else
        msk = imread(msk_filename);
        msk = msk > 0; % in case the mask was produced externally (imageJ)
    end
end