function I = maxProjection(stk)
% MAXPROJECTION   Returns the maximum projection (along first dimension) of an image stack
%
%   @input: stk - an image stack with dimensionality of 3 or higher
%
%   @output: I - a single image which takes the maxmimum value across the
%   first dimension of the array. Has dimensionality one fewer than stk.
   
    if(ndims(stk) < 3)
        error('Image stack has fewer than 3 dimensions. Unable to project to single image.');
    end
    I = squeeze(max(stk,[],1));
end