function lbl = refinedLabels(skel)
% REFINEDLABELS  Improved version of LABELMATRIX
%
%   @input: skel - binary skeleton image
%
%   @output: lbl - labelled skeleton image

    lbl = basicLabels(skel);
    bp_lst = getCoordList(getBranchpoints(skel));
    [bp_len, ~] = size(bp_lst);
    indices = [];
    forgotten = [];
    p = getParams();
    %crop_range = 5; %crops an (2n+1) square around bp to find ends + angles
    
    %We keep the two label and three label cases separate otherwise bugs
    %occur due to partial changes
    for idx = 1:bp_len
        i = bp_lst(idx,2); j = bp_lst(idx,1);
        for a = i-1:i+1
            for b = j-1:j+1
                if(lbl(a,b) > 0)
                    indices = [indices lbl(a,b)];
                end
            end
        end
        indices = unique(indices, 'stable'); %ensure no dupilcates
        if(length(indices) == 2)
           lbl = twoFilamentCase(lbl, indices, bp_lst(idx,:));
        end
        indices = [];
    end
    
    for idx = 1:bp_len
        i = bp_lst(idx,2); j = bp_lst(idx,1);
        for a = i-1:i+1
            for b = j-1:j+1
                if(lbl(a,b) > 0)
                    indices = [indices lbl(a,b)];
                end
            end
        end
        indices = unique(indices, 'stable'); %ensure no dupilcates
        if(length(indices) == 3)
            %all info from st# is lost after this, so it doesn't need to be
            %accurate - left over from old stuff
            [lbl_crop, bp_adj] = cropAndAdjust(lbl, bp_lst(idx,:), p.relabelRange);
            st1 = getBranchAngle(indices(1), indices(2), bp_adj, lbl_crop);
            st2 = getBranchAngle(indices(2), indices(3), bp_adj, lbl_crop);
            st3 = getBranchAngle(indices(3), indices(1), bp_adj, lbl_crop);
            st_new = getMinimumAngle(st1,st2,st3);
            idx1 = st_new.filamentIndices(1); idx2 = st_new.filamentIndices(2);
            lbl = relabel(lbl, idx1, idx2);
            lbl(i,j) = idx1; %insert branchpoint to join up the two segments
            forgotten = [forgotten; idx2];
        end
        indices = [];
    end
    forgotten = sort(forgotten);
    fil_start = max(lbl(:));
    ctr = 1; ctr_max = length(forgotten);
    
    for f = fil_start:-1:1
        if(ctr > ctr_max) %finished all 'forgotten' indices
            break;
        elseif(f < forgotten(ctr)) % No gaps left
            break;
        end
        
        if(ismember(f, lbl)) %filament with number 'f' found
            lbl(lbl == f) = forgotten(ctr); % change label
            ctr = ctr + 1;
        end
    end
    lbl = cleanSkeleton(lbl);
    lbl = fixMissingLabels(lbl); %Final check for missing labels
end

% Crops the whole image to look at range (default 5px) away from branch
% point. This is for measuring filament angles and determining which
% pairing are considered 'straight' and which is the 'branch'. Takes care
% of the shifting branchpoint
function [lbl_crop, bp_adj] = cropAndAdjust(lbl, bp, range)
    x = bp(1); y = bp(2);
    [sizeY,sizeX] = size(lbl);
    [startX,endX] = getRangeLims(x,range,sizeX);
    [startY,endY] = getRangeLims(y,range,sizeY);
    lbl_crop = lbl(startY:endY, startX:endX);
    bp_adj(1) = x - startX + 1; %+1 is because index starts at 1
    bp_adj(2) = y - startY + 1;  
end

function lbl = twoFilamentCase(lbl, indices, bp)
    skel_two = (lbl == indices(1) | lbl == indices(2)); %binary mask of two filaments
    lbl_two = lbl; lbl_two(~skel_two) = 0; %apply mask to return only two labelled filament
    
    lbl_three = basicLabels(rmFourConn(lbl_two,bp) > 0); %remove 4-connected pixels to fragment into 3 filaments
    if(max(lbl_three(:)) ~= 3) %some very exotic shapes can cause problems and are ignored
        return;
    end
    points_lost = getCoordList(xor(lbl_three > 0, skel_two));
    if(size(points_lost,1) ~= 2) %Catches edge cases that should be left alone, such as a loop
        return;
    end
    lbl_three = addPoints(lbl_three,points_lost); % Add back in the points lost in rmFourConn()
    [lbl_three_crop, bp_adj] = cropAndAdjust(lbl_three, bp, 5);
    st1 = getBranchAngle(1,2, bp_adj, lbl_three_crop);
    st2 = getBranchAngle(2,3, bp_adj, lbl_three_crop);
    st3 = getBranchAngle(3,1, bp_adj, lbl_three_crop);
    st_new = getMinimumAngle(st1,st2,st3);
    idx1 = st_new.filamentIndices(1); idx2 = st_new.filamentIndices(2);
    lbl_three = relabel(lbl_three, idx1, idx2);
    lbl_three(bp(2),bp(1)) = idx1; %insert branchpoint to join up the two segments
    branch_label_num = getBranchLabelNumber(idx1,idx2);
    branch_skel = (lbl_three == branch_label_num);
    filament_skel = (lbl_three == idx1);
    lbl(filament_skel) = indices(1); %the number doesn't matter, as long as they are labelled correctly with the two different values from indices
    lbl(branch_skel) = indices(2);
end

% Sets pixels adjcent North, East, South and West of specified point to 0
function lbl = rmFourConn(lbl, p)
    lbl(p(2)+1,p(1)) = 0; %North
    lbl(p(2),p(1)+1) = 0; %East
    lbl(p(2)-1,p(1)) = 0; %South
    lbl(p(2),p(1)-1) = 0; %West
end

function lbl = addPoints(lbl, points)
    len = size(points,1);
    idx_to_add = zeros(len,1);
    for idx = 1:len
        p = points(idx,:); %look around removed points to find the label number
        for i = p(2)-1:p(2)+1
            for j = p(1)-1:p(1)+1
                if(lbl(i,j) ~= 0)
                    idx_to_add(idx) = lbl(i,j);
                end
            end
        end
    end
    for idx = 1:len
        lbl(points(idx,2),points(idx,1)) = idx_to_add(idx); %add back in the points with the correct label
    end
end

%This is a quick hacky way of doing it that is slow to execute - replace
function n = getBranchLabelNumber(idx1, idx2)
    for i = 1:3
        if(idx1 ~= i && idx2 ~= i)
            n = i;
        end
    end
end

% This function removes branchpoints so that filament lengths between
% branchpoints can be determined.
% It removes it twice in the case where the first removal still leaves
% diagonal (8-connectivity) connections.
% Returns a labeled skeleton
function skel_label = basicLabels(im_skel)
    im_skel_seg = im_skel - bwmorph(im_skel, 'branchpoints');
    im_skel_seg = im_skel_seg - bwmorph(im_skel_seg, 'branchpoints'); %ensure that no labels have three ends...
    im_skel_seg = bwmorph(im_skel_seg, 'clean');    %ensure no lone pixels
    %im_skel_seg = im_skel - getBranchpoints(im_skel);
    skel_label = labelmatrix(bwconncomp(im_skel_seg));
end

function l = relabel(l, idx_new, idx_old)
% RELABEL changes one label index to another in labelled image
%
%   @input: l - original labelled image
%           idx_new - new index
%           idx_old - old index
%
%   @output: l - label matrix with updated labelling

    l(l == idx_old) = idx_new;
end

function ba = getBranchAngle(idx1, idx2, bp, im_label)
    im_label(im_label ~= idx1 & im_label ~= idx2) = 0; %remove all other filaments we don't care about
    im_ends = bwmorph((im_label > 0), 'endpoints');

    for i = bp(2)-1:bp(2)+1
        for j = bp(1)-1:bp(1)+1
            im_ends(i,j) = 0; % extra endpoints touching branchpoint, we need to remove them
        end
    end
    ends = getCoordList(im_ends);
    
    
    if(size(ends,1) == 2) %ensure two coordinate pairs
        ba.branchPoint = bp;
        ba.filamentIndices = [idx1,idx2];
        ba.ends = ends;
        ba.angXY = pi - angleMapNew(ends(1,:),ends(2,:),bp);
    else
        ba.angXY = NaN;
    end
end

function ba = getMinimumAngle(ba1, ba2, ba3)
    angles = [ba1.angXY;ba2.angXY;ba3.angXY];
    cor_ang = min(angles);
    if(cor_ang == ba1.angXY)
        ba = ba1;
    elseif(cor_ang == ba2.angXY)
        ba = ba2;
    else
        ba = ba3;
    end
end

function ang = angleMapNew(ep1, ep2, bp)
    fil_len1    = norm(bp-ep1);
    fil_len2    = norm(bp-ep2);
    fil_len_res    = norm(ep2-ep1);
    %fil_len1    = pdist([bp;ep1],'euclidean'); % pdist requires machine learning toolbox
%     fil_len2    = pdist([bp;ep2],'euclidean');
%     fil_len_res = pdist([ep2;ep1],'euclidean');
    cos_rule = (fil_len1.^2 + fil_len2.^2 - fil_len_res.^2) / (2*fil_len1*fil_len2);
    ang = acos(cos_rule);
end

%Make sure each individual filament satisfies the skeleton requirements
function lbl = cleanSkeleton(lbl)
    for i = 1:max(lbl(:))
        curr_fil = lbl==i;
        lbl(xor(curr_fil, bwskel(curr_fil)))=0; %px outside skel requirements are deleted
        if(sum(curr_fil(:)) < 2)
            lbl(curr_fil) = 0;
        end
    end
end

function lbl = fixMissingLabels(lbl)
    used_num = unique(lbl);
    max_num = max(lbl(:));
    unused_num = setdiff(0:max_num,used_num); %start from zero as that is used for background
    len = length(unused_num);
    if(len>0) %only relabel if needed
        for i = 1:length(unused_num)
            high_num = max_num - i + 1;
            if(high_num > unused_num(i))
                lbl(lbl==high_num) = unused_num(i);
            end
        end
    end
end