function ang = angleMap(ep1, ep2, bp)
%ANGLEMAP  Finds the angle between the lines ep1->bp and bp->ep2
%
%   @input: ep1 - Starting point for travel
%           ep2 - End point for travel
%           bp  - Central waypoint between ep1 and ep2
%   
%   @output: ang - The angle of the triangle ep1->bp->ep2
%
% The cosine rule is used to determine the angle
    fil_len1 = norm(bp-ep1);
    fil_len2 = norm(bp-ep2);
    fil_len_res = norm(ep2-ep1);
%     fil_len1    = pdist([bp;ep1],'euclidean');
%     fil_len2    = pdist([bp;ep2],'euclidean');
%     fil_len_res = pdist([ep2;ep1],'euclidean');
    cos_rule = (fil_len1.^2 + fil_len2.^2 - fil_len_res.^2) / (2*fil_len1*fil_len2);
    ang = acos(cos_rule);
end