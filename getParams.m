function p = getParams()
% GETPARAMS  Returns a structure containing all the parameters used in NetEx.
%
%   @output: p - structure containing all used parameters
%
%   Each parameter has a little description explaining their meaning and
%   usage.

    p.z_scale = 5.5; %ratio of pixel size to voxel depth
    p.rb_radius = 15; % 'radius' of rolling ball style BG filtering
    p.curve_range = 6; %px length scale for curvature from a point (so total len is doubled)
    p.prctile_thresh = 80; %percentile for binary image threshold
    p.minBranchLength = 10; %minimum length of branches in a skeleton (to remove 'hairs')
    p.minBlob = 20; %minimum number of connected pixels in binary during skeleton refinement
    p.relabelRange = 5; %distance down each filament to go to find filament angles in relabelling algorithm
    p.widthRange = 2; %px distance used to calc filament direction, so normal can be found for width measurement
    p.maxWidth = 14; %maximum px of width measurements, useful for busy areas (an even int is a good choice)
    p.widthMinBinary = 50; %Minimum number of connected pixels in binary used for width measurements
end
