function I = maskImage(I,msk)
% MASKIMAGE  Sets all pixels outside of mask to zero
%
%   @input: I - image to be masked
%           msk - mask to be applied
%
%   @output: I - masked image
%
%   Handles non-binary mask types (e.g. from ImageJ)

    msk = msk > 0;
    I(~msk) = 0;
end