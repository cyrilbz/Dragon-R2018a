function c = getCurvatureSigned(lbl, range)
% GETCURVATURESIGNED  Calculates the mean curvature for labelled filaments
% by summing the curvature from 3 pixels on a filament separated by a
% designated range. This is repeated for all non-end pixels, accounting for
% curvature pointing inside or outside the shape or path of the filament.
% Absolute values are used for the final result.
%
%   @input: lbl - Labelled skeleton, such as those returned by labelmatrix
%           range - distance between pixels for curvature calculation
%
%   @output: c - Array of mean curvatures, with idx matching filament label


    len = max(lbl(:));
    c = zeros(1,len);
    for l = 1:len
        filament_skel = bwskel(lbl == l);
        inside_bin = createInsideBinary(filament_skel);
        coords = getFilamentOrder(filament_skel);
        if(isnan(coords)) %Very rare case of filaments not working, insert NaN and move on.
            c(l) = NaN;
            continue;
        end
        sum = 0.0; fil_len = size(coords,1);
        if(fil_len > 2)
            for i = 2:fil_len-1 %ignore endpoints for curvature
                point = coords(i,:);
                [idxMin,idxMax] = getRangeLims(i,range,fil_len); %coordinate idx of ends for curvature
                e1 = coords(idxMin,:); e2 = coords(idxMax,:);
                curve_sign = 1 ;%getCurveSign(e1,e2,point,inside_bin);
                curve = curve_sign * measureCurve(e1,e2,point);
                sum = sum + curve;
            end
            c(l) = abs(sum/(fil_len-2));
        else
            c(l) = 0;
        end
    end
end

function [e1,e2] = getEnds(p, skel, range)
    [sizeY,sizeX] = size(skel);
    [startX,endX] = getRangeLims(p(1), range, sizeX);
    [startY,endY] = getRangeLims(p(2), range, sizeY);

    %sub_skel = skel(startY:endY, startX:endX);
    skel_masked = zeros(size(skel));
    skel_masked(startY:endY, startX:endX) = skel(startY:endY, startX:endX);
    ends = getCoordList(bwmorph(skel_masked,'endpoints'));
    e1 = ends(1,:); e2 = ends(2,:);
end

function sign = getCurveSign(e1,e2,p,bin)
    x_ends = [e1(1); e2(1)];
    y_ends = [e1(2); e2(2)];
    p_x = p(1);
    if(x_ends(1) == x_ends(2))
        p_y_interp = round(mean(y_ends(:)));
    else
        p_y_interp = round(interp1(x_ends,y_ends,p_x));
        if(isnan(p_y_interp))
            p_y_interp = p(2);
            if(y_ends(1) == y_ends(2))
                p_x = round(mean(x_ends(:)));
            else
                p_x = round(interp1(y_ends,x_ends,p_y_interp));
            end
        end   
    end
    if(bin(p_y_interp,p_x))
        sign = 1;
    else
        sign = -1;
    end
end


function bin = createInsideBinary(s)
    ends = getCoordList(bwmorph(s,'endpoints'));
    if(size(ends,1) < 2) %less than two end points = closed loop
        bin = imfill(s,'holes');
        return;
    end
    bin = s;
    for i = 1:size(ends,1)
        x = ends(i,1); y = ends(i,2);
        dx = 0; dy = 0;
        for a = -1:1
            for b = -1:1
                if(s(y+b,x+a))
                    if(a ~= 0 || b ~= 0)
                        dx = -a; dy = -b; %switch sign to go from adj to end
                    end
                end
            end
        end
        [sizeY,sizeX] = size(s);
        while(x > 1 && x < sizeX && y > 1 && y < sizeY)
            x = x + dx; y = y + dy;
            bin(y,x) = 1;
        end
    end
    if(size(getCoordList(bwmorph(bin,'endpoints')),1) == 1) %closed loop, can fill
        bin = imfill(bin, 'holes');
    else
        mean_point = pointsToFill(bin);
        bin = imfill(bin,flip(mean_point)); %flip as coords are [x,y] but idx is [y,x]
    end
end
function p = pointsToFill(bin)
    new_ends = getCoordList(bwmorph(bin,'endpoints'));
    p = round(mean(new_ends));
    if(bin(p(2),p(1)))
        p = round(mean(getCoordList(bin)));
    end
end