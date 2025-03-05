function I8 = to8Bit(I)
% TO8BIT  Converts an image to 8-bit, scaling values to fill 8-bit range
%
%   @input: I - an image to be converted to 8-bit
%
%   @output: I8 - the 8-bit version of the input, I

    [sizeY,sizeX] = size(I);
    max8bit = 255;
    min_val = min(I(:));
    im_shift = I - min_val;
    max_val = max(im_shift(:));
    scale_fctr = double(max8bit)/double(max_val);
    
    I8 = im_shift;
    for i = 1:sizeY
        for j = 1:sizeX
            I8(i,j) = round(scale_fctr * im_shift(i,j));
        end
    end
    I8 = uint8(I8); %Ensure Matlab knows the data type is now uint8
end