function [idx_start,idx_end] = getRangeLims(idx,range,max)
% GETRANGELIMS  returns the start and end indicies for a search of length
% range in both directions from central pixel idx. Accounts for minumum and
% maximum values
%
%   @input: idx   - central pixel to start the ranged search from
%           range - distance in pixels from idx to look at. Total distance
%                   up to (2*range) + 1
%           max   - maximum value the range can look at, such as the edge
%                   of the image
%
%   @output: idx_start - the minimum value of the search, >= (idx - range)
%            idx_end   - the maximum value of the search, <= (idx + range)

    idx_start = idx - range; idx_end = idx + range;
    if(idx_start < 1)
        idx_start = 1;
    end
    if(idx_end > max)
        idx_end = max;
    end
end