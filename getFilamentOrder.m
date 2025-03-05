function coords = getFilamentOrder(s)
% GETFILAMENTORDER  Returns co-ordinates of true points from a skeletonised
% image, in order, starting from the left-most endpoint (or uppermost in
% the case of a tie).
%
%   @input: s - skeletonised image. If the filaments are >1px wide or there
%   are more than 1 filament, this function will not work
%
%   @output: coords - An Nx2 array of co-ordinates in the skeleton

    N = sum(s(:));
    coords = zeros(N,2);
    ends = getCoordList(bwmorph(s,'endpoints'));
    if(size(ends,1) == 0)
        c = getCoordList(s);
        ends = c(1,:); %grab a random point to start on in case of a loop
    end
    coords(1,:) = ends(1,:);
    x0 = ends(1,1); y0 = ends(1,2);
    ctr = 2;
    while(ctr <= N)
        skel_pts = [];
        for y = y0-1:y0+1
            for x = x0-1:x0+1
                if(s(y,x))
                    skel_pts = [skel_pts;x,y]; %grab all points connected to (x0,y0)
                end
            end
        end
        new_pts = setdiff(skel_pts,coords,'stable','rows'); %skeleton not already accounted for
        if(size(new_pts,1) == 2) %Need to pick 4-connected point first, then 8-connected
            if(new_pts(1,1) == x0 || new_pts(1,2) == y0) %4-connected if one of these are true
                coords(ctr,:) = new_pts(1,:);
                coords(ctr+1,:) = new_pts(2,:);
                x0 = new_pts(2,1); y0 = new_pts(2,2);
            else
                coords(ctr,:) = new_pts(2,:);
                coords(ctr+1,:) = new_pts(1,:);
                x0 = new_pts(1,1); y0 = new_pts(1,2);
            end
            ctr = ctr + 2;
        elseif(size(new_pts,1) == 1)
            coords(ctr,:) = new_pts(1,:);
            x0 = new_pts(1,1); y0 = new_pts(1,2);
            ctr = ctr + 1;
        else
            warning("Unable to find new pixel despite not having found them all. Likely a loop within a loop.");
            coords = NaN;
            break; %some skeletonisation issue
        end
    end
end