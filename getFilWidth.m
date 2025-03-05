function fw = getFilWidth(im,skel_lbl,msk)
% GETFILWIDTH  Estimates the thickness of filaments from a labelled skeleton
%
%   @input: im       - Raw image used to measure filament thickness
%           skel_lbl - Labelled skeleton derived from im
%           msk      - Mask around a single cell in the image
%
%   @output: fw - A list of pixel widths, with index matching the labels of the skeleon

    p = getParams();
    %range = 2; %distance from centre to find end point to look at local direction. Width is perpendicular to this direction
    wb = widthBinary(im,msk);
    fils = max(skel_lbl(:));
    fw = zeros(fils,1);
    [sizeY, sizeX] = size(skel_lbl);
    
    for i = 1:fils %loop over different filaments in region
        filIm = (skel_lbl == i);
        [coY,coX] = find(filIm); % gets x and y coordinates of i'th filament
        sumWid = 0; len = length(coY);
        
        prev = 0;
        for j = 1:len %loop over the pixels in a single filament
            cenY = coY(j); cenX = coX(j); % co-ords of pixel to look at
            [startX,endX] = getRangeLims(cenX,p.widthRange,sizeX); %local region to determine filament direction
            [startY,endY] = getRangeLims(cenY,p.widthRange,sizeY);
            widX = endX - startX; widY = endY - startY; 
            seg = imcrop(filIm, [startX startY widX widY]); %image of local region
            [endY,endX] = find(bwmorph(seg, 'endpoints')); %two endpoints for local region
            dy = endY(2) - endY(1); dx = endX(2) - endX(1);
            ang = atan2(dy,dx);
            curr = getWidth(wb,cenX,cenY,ang);
            if(prev ~= 0) %Don't allow width to change by more than two pixels for next px in filament
                if(curr > prev + 2)
                    curr = prev + 2;
                elseif(curr < prev - 2)
                    curr = prev - 2;
                end
            end
            sumWid = sumWid + curr;
            prev = curr;
        end
        fw(i) = double(sumWid)/double(len);
    end
end
function wb = widthBinary(im, msk)
    p = getParams();
    im(~msk) = 0; % mask image
    fib = fibermetric(im);
    m = meanMasked(fib,msk);
    wb = bwareaopen(fib > m,p.widthMinBinary);
end
function w = getWidth(im, cenX, cenY, dir)  %dir is direction of filament
    [sizeY, sizeX] = size(im);
    pi2 = pi/2.0;
    p = getParams();
    lim = p.maxWidth/2;
    dir1 = dir + pi2; dir2 = dir - pi2; %rotate dir by 90deg for normals
    
    for i = 1:lim
        x = cenX + round(i*cos(dir1)); y = cenY + round(i*sin(dir1));
        if(x < 1 || x > sizeX || y < 1 || y > sizeY)
            break;
        elseif(im(y,x) == 0)
            break;
        end
        
    end
    for j = 1:lim
        x = cenX + round(j*cos(dir2)); y = cenY + round(j*sin(dir2));
        if(x < 1 || x > sizeX || y < 1 || y > sizeY)
            break;
        elseif(im(y,x) == 0)
            break;
        end
    end
    w = i+j-1; % central pixel is counted twice, -1 accounts for this
end