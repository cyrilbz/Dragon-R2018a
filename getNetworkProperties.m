function np = getNetworkProperties(im_raw, im_mask)
% GETNETWORKPROPERTIES  This function extracts the properties of an actin network from an image
%
%   @input: im_raw  - a 2D image or projection of the actin network
%           im_mask - a binary mask around the cell of interest
%
%   @output: np - A range of network properties in a structure

    p = getParams();
    
    np.orientation = getOrientation(im_mask);
    im_raw_rt = imrotate(im_raw,rad2deg(-np.orientation));
    im_mask2 = imrotate(im_mask, rad2deg(-np.orientation));
    np.imRaw = maskImage(im_raw_rt,im_mask2);
    np.imMask = im_mask2;
    ini_skel = multiOrientSkel(im_raw,im_mask,p.prctile_thresh);
    bin = getBinary(im_raw, im_mask, p.prctile_thresh);
    np.imSkel = refineSkeleton(ini_skel);
    disp('Skeleton done')
    np.cellSize = sum(im_mask2(:));
    np.skelDensity = sum(np.imSkel(:))/double(np.cellSize);
    np.actDensity = sum(bin(:))/double(np.cellSize) ;
    branchpoints = bwmorph(np.imSkel, 'branchpoints');
    np.cellBPDensity = numel(branchpoints)/ double(np.cellSize);
    np.Structures = structureSize(np.imSkel);
    np.Structures.freq = double(length(np.Structures.sizes)) / double(np.cellSize);
    np.skelLabel = refinedLabels(np.imSkel);
    disp('Labels done - starting measurements')
    np.filWidth = getFilWidth(np.imRaw, np.skelLabel, im_mask2);
    disp('Width : done')
    [np.deviation, np.newcurv] = avgDeviation(np.skelLabel);
    disp('Deviation : done')
    %np.curvature = getCurvature(np.skelLabel, p.curve_range);
    np.curvatureSigned = getCurvatureSigned(np.skelLabel, p.curve_range);
    disp('Curvatures : done')
    np.filLenXY = getFilamentLengths(np.skelLabel);
    np.avgLenMapXY = avgFilamentLengthMap(np.skelLabel, np.filLenXY);
    np.branchDen = branchDensity(np.imSkel);
    np.actinDen = actinDensity(np.imSkel);
    np.filamentAng = filamentAngles(np.skelLabel);
    np.branchAng = RelabelAngles(np.skelLabel);
    np.branchRatio = double(size(np.branchAng))/double(size(np.filamentAng));
end

function s = refineSkeleton(s)
% REFINESKELETON  Dilates and erodes skeleton to remove small loops & hairs
%
%   @input: s - binary skeleton
%
%   @output: s - refined binary skeleton
    
    p = getParams();
    s = bwmorph(s, 'close');
    s = bwskel(s, 'MinBranchLength', p.minBranchLength);
    s = bwareaopen(s,p.minBlob);
end

function SS = structureSize(imSkel)
% STRUCTURESIZE  determines the mean and stdev of the size of connected 
% actin structures
%
%   @input: imSkel - a binary skeleton
%
%   @output SS - structure containing the mean and stdev of structure sizes
%
%   This uses connectivity of the skeleton to determine distinct structures

    l = labelmatrix(bwconncomp(imSkel));
    N = max(l(:));
    SS.sizes = zeros(N, 1);
    for i = 1:N
        buff = l==i;
        SS.sizes(i) = sum(buff(:));
    end
    
    SS.mean = mean(SS.sizes);
    SS.stdev = std(SS.sizes);
end

%This function measures the average distance between the filament and a
%straight line which goes between the two endpoints.
function [d,k] = avgDeviation(skel_label)
    no_of_filaments = max(skel_label(:));
    d = zeros(1, no_of_filaments);
    k = zeros(1, no_of_filaments); 
    for f = 1:no_of_filaments
        skel = (skel_label == f);
        coords = getCoordList(skel);
        ends = getCoordList(bwmorph(skel, 'endpoints'));
        if(size(ends,1) == 2)
            sum_d = 0.0; % for deviation
            len = size(coords,1);
            h = zeros(1,len); % to save the distance
            kappa = zeros(1,len); % to save the curvature
            for j = 1:len % compute the deviation
                h(j) = distToLine(ends, coords(j,:));
                sum_d = sum_d + h(j); % sum all deviations
            end
            %  ! the part below could be turned into a matrix operation !
            for j = 2:len-1 % compute the curvature using height-function
                hx = 1/2*(h(j+1)-h(j-1)) ; % centered 1st derivative
                hxx = h(j+1)+h(j-1)-2*h(j) ; % centered 2nd derivative
                kappa(j) = hxx/(1+hx^2)^3/2 ; % curvature (signed)
            end
            % compute curvature at boundaries 
            % (extrapole at constant values - 1st order?)
            hx = 1/2*(h(2)-h(1)) ; hxx =  h(2)-h(1) ;
            kappa(1) = hxx/(1+hx^2)^3/2 ;
            hx = 1/2*(h(len)-h(len-1)) ; hxx =  h(len-1)-h(len) ;
            kappa(len) = hxx/(1+hx^2)^3/2 ;
            
            % compute mean values
            d(f) = sum_d / double(len); % mean deviation
            k(f) = sum(abs(kappa))/double(len) ; % mean absolute curvature
        end
    end
end
%{
function r = distToLine(ends, p) %could be replaced by pdist?
    dy = ends(1,2) - ends(2,2); dx = ends(1,1) - ends(2,1);
    if(dx == 0)
        r = abs(p(1) - ends(1,1));
    elseif(dy == 0)
        r = abs(p(2) - ends(1,2));
    else
        line_len = sqrt(dy.^2 + dx.^2);
        e1 = ends(1,:); e2 = ends(2,:);
        r = abs((dy*p(1)) - (dx*p(2)) + (e2(1)*e1(2)) - (e2(2)*e1(2)))/line_len;
    end
end
%}
function r = distToLine(ends,p)
    e1 = ends(1,:); e2 = ends(2,:);
    A = triArea(e1,e2,p);
    %line_len = pdist2(e1,e2,'euclidean');
    line_len = norm(e1-e2);
    r = (2*A)/line_len;
end

% This function calculates the average Menger Curvature for each filament
% It returns an array of curvatures for each labelled filament
function c = getCurvature(skel_label, range)
    len = max(skel_label(:));
    c = zeros(1,len);
    for l = 1:len
        filament_skel = (skel_label == l);
        ends = bwmorph(filament_skel,'endpoints');
        centre_fil = filament_skel - ends;
        centre_coords = getCoordList(centre_fil);
        sum = 0.0; pxCentreLen = size(centre_coords,1);
        if(pxCentreLen > 0)
            for i = 1:pxCentreLen
                point = centre_coords(i,:);
                [e1,e2] = getEnds(point, filament_skel, range);
                sum = sum + measureCurve(e1,e2,centre_coords(i,:));
            end
            c(l) = sum/pxCentreLen;
        else
            c(l) = 0;
        end
    end
end

function [e1,e2] = getEnds(p, skel, range)
    [sizeY,sizeX] = size(skel);
    [startX,endX] = getRangeLims(p(1), range, sizeX);
    [startY,endY] = getRangeLims(p(2), range, sizeY);

    sub_skel = skel(startY:endY, startX:endX);
    ends = getCoordList(bwmorph(sub_skel,'endpoints'));
    e1 = ends(1,:); e2 = ends(2,:);
end

% This function calculates the average length of filaments in proximity of
% every pixel (range is arbitrary)
% Returns a matrix of the average filament lengths
function avg_len_map = avgFilamentLengthMap(skel_lbl, fil_len)
    [sizeY, sizeX] = size(skel_lbl);
    avg_len_map = zeros(size(skel_lbl));
    range = 10;
    coords = getCoordList(skel_lbl > 0);
    no_of_coords = size(coords,1);
    for c = 1:no_of_coords
        x = coords(c,1); y = coords(c,2);
        [startY,endY] = getRangeLims(y,range,sizeY);
        [startX,endX] = getRangeLims(x,range,sizeX);
        no_of_pixels = (endY - startY + 1)*(endX - startX + 1);
        density_increase = fil_len(skel_lbl(y,x))/no_of_pixels;
        for i = startY:endY
            for j = startX:endX
                avg_len_map(i,j) = avg_len_map(i,j) + density_increase;
            end
        end
    end
    avg_len_map(avg_len_map == 0) = NaN;
end

% This function calculates the density of branch points in proximity of
% each pixel (range is arbitrary)
% Returns a matrix of these densities for each pixel
function branch_den = branchDensity(im_skel)
    [sizeY, sizeX] = size(im_skel);
    branch_den = zeros(size(im_skel));
    im_branch = getBranchpoints(im_skel);
    range = 20;
    for i = 1:sizeY
        for j = 1:sizeX
            if(im_branch(i,j))
                [startY,endY] = getRangeLims(i,range,sizeY);
                [startX,endX] = getRangeLims(j,range,sizeX);
                no_of_pixels = (endY - startY + 1)*(endX - startX + 1);
                density_increase = 1/no_of_pixels;
                for y = startY:endY
                    for x = startX:endX
                        branch_den(y,x) = branch_den(y,x) + density_increase;
                    end
                end
            end
        end
    end
    branch_den(branch_den == 0) = NaN; %set NaN where no branchpoints are so that its distinguished from low values
end

% This function looks at the density of pixels in the skeleton in proximity
% of each point (range is arbitrary)
% Returns a matrix of these densities for each pixel
function actin_den = actinDensity(im_skel)
    [sizeY,sizeX] = size(im_skel);
    actin_den = zeros(size(im_skel),'double');
    range = 15;
    for i = 1:sizeY
        for j = 1:sizeX
            if(im_skel(i,j))
                [startY,endY] = getRangeLims(i,range,sizeY);
                [startX,endX] = getRangeLims(j,range,sizeX);
                no_of_pixels = (endY - startY + 1)*(endX - startX + 1);
                density_increase = 1/no_of_pixels;
                for y = startY:endY
                    for x = startX:endX
                        actin_den(y,x) = actin_den(y,x) + density_increase;
                    end
                end
            end
        end
    end
    actin_den(actin_den == 0) = NaN; %set NaN where no actin filaments are so that its distinguished from low values
end

% Ends are stored as [x,y ; x,y]
% This function gets the angle of each branch realtive to the +ve x-axis
function filament_ang = filamentAngles(skel_label)
    % generate array of structures
    no_of_filaments = max(skel_label(:));
    fa.ends = [0,0 ; 0,0]; fa.angXY = 0.0; fa.angZ = 0.0;
    filament_ang = repmat(fa,1,no_of_filaments);
    
    % grab endpoints of each filament
    im_skel = skel_label > 0;
    ep = bwmorph(im_skel, 'endpoints');
    [sizeY,sizeX] = size(skel_label);
    for i = 1:sizeY
        for j = 1:sizeX
            if(ep(i,j) == 1)
                label = skel_label(i,j);
                if(label == 0)
                    error("Error when finding endpoints");
                end
                if(filament_ang(label).ends(1,1) == 0)
                    filament_ang(label).ends(1,1) = j;
                    filament_ang(label).ends(1,2) = i;
                elseif(filament_ang(label).ends(2,1) == 0)
                    filament_ang(label).ends(2,1) = j;
                    filament_ang(label).ends(2,2) = i;
                else
                    error("Filament has more than two ends!");
                end
            end
        end
    end
    
    % Determine angle relative to cell axis
    for i = 1:no_of_filaments
        dy = filament_ang(i).ends(2,2) - filament_ang(i).ends(1,2);
        dx = filament_ang(i).ends(2,1) - filament_ang(i).ends(1,1);
        filAng = atan2(dy,dx);
        if(filAng > pi/2.0)
            filAng = filAng - pi;
        elseif(filAng < -pi/2.0)
            filAng = filAng + pi;
        end
        filament_ang(i).angXY = filAng;
    end
end

function ba_lst = RelabelAngles(lbl)
    p = getParams();
    bp_lst = getCoordList(bwmorph(lbl > 0, 'branchpoints'));
    len = size(bp_lst,1);
    ba_lst = [];
    
    [sizeY,sizeX] = size(lbl);
    for bp_num = 1:len
        cur_bp = bp_lst(bp_num,:);
        [startX,endX] = getRangeLims(cur_bp(1),p.relabelRange,sizeX); %don't go over edges of image in range search
        [startY,endY] = getRangeLims(cur_bp(2),p.relabelRange,sizeY);
        width = endX-startX; height = endY-startY;
        im = imcrop(lbl, [startX startY width height]);
        temp_bpX = cur_bp(1) - startX + 1; %offset by 1 because matlab indexing starts at 1
        temp_bpY = cur_bp(2) - startY + 1;
        ba_lst = [ba_lst getBA(im, temp_bpX, temp_bpY)];
        %adjust bp and ends outside so that getBA() can be used on non-cropped images
        ba_lst(bp_num).branchPoint = cur_bp;
        ba_lst(bp_num).ends = adjustEnds(ba_lst(bp_num).ends, startX, startY);
    end
end
function ba = getBA(lbl, bpX, bpY) 
    indices = [];
    for i = bpY-1:bpY+1
        for j = bpX-1:bpX+1
            if(lbl(i,j) > 0)
                indices = [indices; lbl(i,j)];
            end
        end
    end
    indices = unique(indices,'stable');
    if(length(indices) == 1) %We have a closed loop with a line coming off (like a tennis racquet)
        ba = createEmptyBranchAngStruct();
        ba.branchPoint = [bpX bpY];
        ba.ends = [getCoordList(bwmorph(lbl==indices,'endpoints')); NaN, NaN];
        ba.angXY = NaN;
        ba.angZ = NaN;
        ba.filamentIndices = [indices, NaN];
        return;
    end
    lbl(lbl ~= indices(1) & lbl ~= indices(2)) = 0; %remove all other filaments we don't care about
    
    ep_im = bwmorph(lbl == indices(1), 'endpoints') + bwmorph(lbl == indices(2), 'endpoints'); %find ep of separate filaments and add -> fixes issues with looped branches
    if(sum(ep_im(:)) > 3)
        ep_im(bpY-1:2:bpY+1,bpX) = 0; ep_im(bpY,bpX-1:2:bpX+1) = 0; %remove 4-connected endpoints
    end
    if(sum(ep_im(:)) > 3)
        ep_im(bpY-1:2:bpY+1,bpX-1) = 0; ep_im(bpY,bpX-1:2:bpX+1) = 0; %remove remaining 4 endpoints diagonally connected to branchpoint
    end
    ep_lst = getCoordList(ep_im);
    
    ba = determineBA(ep_lst, bpX, bpY);
    ba.filamentIndices(1) = indices(1); ba.filamentIndices(2) = indices(2);
    ba.branchPoint(1) = bpX; ba.branchPoint(2) = bpY;
end
function ba = determineBA(ep_lst, bpX, bpY)
    ba1 = createEmptyBranchAngStruct();
    ba2 = createEmptyBranchAngStruct();
    ba3 = createEmptyBranchAngStruct();
    
    ba1.ends(1,:) = ep_lst(1,:); ba1.ends(2,:) = ep_lst(2,:);
    ba1.angXY = pi - angleMap(ba1.ends(1,:), ba1.ends(2,:), [bpX bpY]);

    ba2.ends(1,:) = ep_lst(2,:); ba2.ends(2,:) = ep_lst(3,:);
    ba2.angXY = pi - angleMap(ba2.ends(1,:), ba2.ends(2,:), [bpX bpY]); %(pi - ang) gives the angle of deviation for anything travelling along cable
    
    ba3.ends(1,:) = ep_lst(3,:); ba3.ends(2,:) = ep_lst(1,:);
    ba3.angXY = pi - angleMap(ba3.ends(1,:), ba3.ends(2,:), [bpX bpY]);

    ba = getCorrectBranchAngle(ba1, ba2, ba3);
end
%branchpoint,ends,anyXY,angZ,filamentIndices
function ba = createEmptyBranchAngStruct()
    ba.branchPoint = [0 0];
    ba.ends = [0 0; 0 0];
    ba.angXY = -1.0;
    ba.angZ = -1.0;
    ba.filamentIndices = [0 0];
end
function ends = adjustEnds(ends,startX,startY)
    ends(1,1) = ends(1,1) + startX - 1; % -1 is because startX was shifted to 1, not 0, so this is accounted for
    ends(1,2) = ends(1,2) + startY - 1;
    ends(2,1) = ends(2,1) + startX - 1;
    ends(2,2) = ends(2,2) + startY - 1;
end

function ang = angleMap(ep1, ep2, bp)
    fil_len1 = norm(bp-ep1);
    fil_len2 = norm(bp-ep2);
    fil_len_res = norm(ep2-ep1);
%     fil_len1    = pdist([bp;ep1],'euclidean');
%     fil_len2    = pdist([bp;ep2],'euclidean');
%     fil_len_res = pdist([ep2;ep1],'euclidean');
    cos_rule = (fil_len1.^2 + fil_len2.^2 - fil_len_res.^2) / (2*fil_len1*fil_len2);
    ang = acos(cos_rule);
end


function ba = getCorrectBranchAngle(ba1, ba2, ba3)
    angles = [ba1.angXY;ba2.angXY;ba3.angXY];
    cor_ang = median(angles);
    if(cor_ang == ba1.angXY)
        ba = ba1;
    elseif(cor_ang == ba2.angXY)
        ba = ba2;
    else
        ba = ba3;
    end
end