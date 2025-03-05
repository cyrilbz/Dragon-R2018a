function lm = localMean(I,range)
% LOCALMEAN  Returns array showing mean in (-range,+range) for every pixel
%
%   @input: I - greyscale image from which to find the local mean
%           range - distance (in pixels, from central pixel) to average
%
%   @output: lm - array of means for every pixel. Same size as I

    [sizeY,sizeX] = size(I);
    lm = zeros(sizeY,sizeX,'double');
    for i = 1:sizeY
        [yStart,yEnd] = getRangeLims(i,range,sizeY);
        for j = 1:sizeX
            [xStart,xEnd] = getRangeLims(j,range,sizeX);
            sum = double(0); pixels = uint32(0);
            for y = yStart:yEnd
                for x = xStart:xEnd
                    val = double(I(y,x));
                    sum = sum + val;
                    pixels = pixels + 1;
                end
            end
            lm(i,j) = sum/double(pixels);
        end
    end
end