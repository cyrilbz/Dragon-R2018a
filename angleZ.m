function net_props = angleZ(net_props, im_stack, z_scale)
% ANGLEZ  Extends the properties from getNetworkProperties() into three dimensions
%
%   @input: net_props - the struct returned from 'getNetworkProperties()'
%           im_stack  - the rotated stack of images forming the 3D network
%           z_scale   - the ratio of pixel height/width to voxel depth
%
%   @output: netProps - extended properties of the extracted actin network

    net_props.filamentAng = filamentAngZ(net_props.filamentAng, im_stack, z_scale);
    [net_props.avgFilXY,net_props.avgFilZ] = meanAngle(net_props.filamentAng);
    [net_props.stdevFilXY,net_props.stdevFilZ] = stdevAngle(net_props.filamentAng,net_props.avgFilXY,net_props.avgFilZ);
    net_props.branchAng = getAngZ(net_props.branchAng, im_stack, z_scale);
    [net_props.avgBranchXY,net_props.avgBranchZ] = meanAngle(net_props.branchAng);
    [net_props.stdevBranchXY,net_props.stdevBranchZ] = stdevAngle(net_props.branchAng,net_props.avgBranchXY,net_props.avgBranchZ);
    [net_props.avgAngMapXY, net_props.avgAngMapZ] = avgAngleMap(net_props.imSkel, net_props.branchAng);
    net_props.filLenXYZ = getFilamentLen3D(net_props.filLenXY, net_props.filamentAng);
    net_props.avgLen = mean(net_props.filLenXYZ(:), 'omitnan');
    net_props.stdevLen = std(net_props.filLenXYZ(:), 'omitnan');
    net_props.avgLenMap = avgFilamentLengthMap(net_props.skelLabel, net_props.filLenXYZ);

    net_props.z_scale = z_scale;
    net_props = volumeDensities(net_props, im_stack);
end

function fil_ang_z = filamentAngZ(filament_ang, im_stack, z_scale)
    for i = 1:length(filament_ang)
        point1 = [filament_ang(i).ends(1,1), filament_ang(i).ends(1,2)];
        point2 = [filament_ang(i).ends(2,1), filament_ang(i).ends(2,2)];
        plane1 = grabPlane(point1, im_stack);
        plane2 = grabPlane(point2, im_stack);
        filament_ang(i).angZ = depthAng(point1, point2, plane1, plane2, z_scale);
    end
    fil_ang_z = filament_ang;
end
function plane = grabPlane(point, im_stack)
    [no_of_planes, maxI, maxJ] = size(im_stack);
    avg_in = zeros(no_of_planes,1);
    x = point(1); y = point(2);
    for p = 1:no_of_planes
        avg = 0; counter = 0;
        for i = y-1:y+1
            for j = x-1:x+1
                if(i > 0 && j > 0 && i < maxI && j < maxJ)
                    avg = avg + im_stack(p,i,j);
                    counter = counter + 1;
                end
            end
        end
        avg_in(p) = avg/counter;
    end
    
    max = 0.0; plane = 0;
    for p = 1:no_of_planes
        if(avg_in(p) > max)
            max = avg_in(p); plane = p;
        end
    end
end
function ang_z = depthAng(bp, ep, plane_bp, plane_ep, z_scale)
    plane_dist = norm(ep-bp);
    %plane_dist = pdist([ep;bp], 'Euclidean');
    height_dist = z_scale * (plane_ep - plane_bp);
    ang_z = atan2(height_dist, plane_dist);
end

function [meanXY,meanZ] = meanAngle(filAng)
    sumXY = 0.0; sumZ = 0.0; N = length(filAng);
    for i = 1:N
        sumXY = sumXY + filAng(i).angXY;
        sumZ  = sumZ  + filAng(i).angZ;
    end
    meanXY = sumXY/N;
    meanZ  = sumZ/N;
end
function [stdevXY,stdevZ] = stdevAngle(filAng,meanXY,meanZ)
    varXY = 0.0; varZ = 0.0; N = length(filAng);
    for i = 1:N
        varXY = varXY + (filAng(i).angXY - meanXY)^2;
        varZ  = varZ  + (filAng(i).angZ  - meanZ)^2;
    end
    stdevXY = sqrt(varXY/N);
    stdevZ  = sqrt(varZ/N);
end

function branch_ang = getAngZ(branch_ang, im_stack, z_scale)
    for i = 1:length(branch_ang)
        ep1 = branch_ang(i).ends(1,:); ep2 = branch_ang(i).ends(2,:);
        bp = branch_ang(i).branchPoint;
        plane1 = getZDist(ep1,im_stack,z_scale); plane2 = getZDist(ep2,im_stack,z_scale);
        bp_plane = getZDist(bp,im_stack,z_scale);
        end1 = [ep1,plane1]; end2 = [ep2 ,plane2];
        branchpoint = [bp,bp_plane];
        branch_ang(i).branchPoint = [bp,bp_plane];
        branch_ang(i).ends = [ep1,plane1 ; ep2,plane2];
        branch_ang(i).angZ = pi - angleMap(end1,end2,branchpoint);
    end
end
function z = getZDist(point,im_stack,z_scale)
    plane_num = grabPlane(point,im_stack);
    z = z_scale * (plane_num - 1);    %Plane number starts from 1, shift to zero before scaling
end

function [avg_ang_map_xy,avg_ang_map_z] = avgAngleMap(im_skel,branch_ang)
    im_branch = bwmorph(im_skel, 'branchpoints');
    im_branch2 = bwmorph(im_skel - im_branch, 'branchpoints');
    im_branch = im_branch + im_branch2;
    
    avg_ang_map_xy = zeros(size(im_skel), 'double');
    avg_ang_map_z  = zeros(size(im_skel), 'double');

    [sizeY,sizeX] = size(im_skel);
    range = 15;
    for i = 1:sizeY
        for j = 1:sizeX
            ang_sum_xy = 0.0; ang_sum_z = 0.0; no_of_angs = 0;
            [startY,endY] = getRangeLims(i,range,sizeY);
            [startX,endX] = getRangeLims(j,range,sizeX);
            for y = startY:endY
                for x = startX:endX
                    if(im_branch(y,x) == 1)
                        [ang_xy, ang_z] = grabAng([x,y], branch_ang);
                        if(ang_xy >= 0 && ang_z >= 0)
                            ang_sum_xy = ang_sum_xy + ang_xy;
                            ang_sum_z  = ang_sum_z  + ang_z;
                            no_of_angs = no_of_angs + 1;
                        end
                    end
                end
            end
            % Allow NaNs to differentiate small angles and no angles present
            avg_ang_map_xy(i,j) = ang_sum_xy/no_of_angs;
            avg_ang_map_z(i,j)  = ang_sum_z/no_of_angs;
        end
    end
end
function [ang_xy,ang_z] = grabAng(bp,branch_ang)
    x = bp(1); y = bp(2);
    no_of_bp = length(branch_ang);
    idx = 1;
    looking_for_angs = true;
    ang_xy = -1.0; ang_z = -1.0; %default values clearly incorrect for error checking
    while(looking_for_angs && idx <= no_of_bp)
        test_x = branch_ang(idx).branchPoint(1);
        if(x == test_x)
            testY = branch_ang(idx).branchPoint(2);
            if(y == testY)
                ang_xy = branch_ang(idx).angXY;
                ang_z  = branch_ang(idx).angZ;
                looking_for_angs = false;
            end
        end
        idx = idx + 1;
    end
end

% This function takes the 2D projection of length and calculates it in 3D
function fil_len_xyz = getFilamentLen3D(fil_len_xy, fil_ang)
    if(length(fil_len_xy) ~= length(fil_ang))
        error("Inconsistent number of actin filaments.");
    end
    fil_len_xyz = zeros(size(fil_len_xy));
    for i = 1:length(fil_len_xyz)
        fil_len_xyz(i) = fil_len_xy(i)/cos(fil_ang(i).angZ);
    end
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