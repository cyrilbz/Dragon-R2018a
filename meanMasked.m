function m = meanMasked(I, msk)
% MEANMASKED  Calculates the mean pixel value in desired region
%
%   @input: I - Image from which the mean will be found
%           msk - Region to calculate the mean in
%
%   @output: m - Mean value of I in region msk

    I(~msk) = NaN;
    %m = mean(I, 'all', 'omitnan'); %requires new versions of Matlab
    m = mean(I(~isnan(I))); % R2018a: Remove NaNs and then calculate the mean
end