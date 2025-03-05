function mdb = meanDistToBranch(branchAng, filLen)
% MEANDISTTOBRANCH  Calcluates the average distance between branchpoints
%
%   @input: branchAng - the branchAng struct(s) obtained from DRAGoN
%           filLen    - the filLen array from DRAGoN. Can use the 2D or 3D
%           versions (filLenXY and filLenXYZ respectively)
%
%   @output: mdb - the mean distance between branchpoints on filaments that
%   branch
%
%   Please note that only the filaments which contain a branch point are
%   included in the average. To get the branch point density relative to
%   the total skeleton length, simply use the following in
%   getNetworkProperties.m:
%       mdb = sum(np.imSkel, 'all')/length(np.branchAng);
%
%   To call this function in getNetworkProperties.m, add the following:
%       np.mdb = meanDistToBranch(np.branchAng, np.FilLenXY);

    bps = length(branchAng);
    branchFils = [];
    for idx = 1:bps
        branchFils = [branchFils branchAng(idx).filamentIndices]; %collect list of filaments that branch
    end
    branchFils = unique(branchFils);
    totalLen = sum(filLen(branchFils), "all");
    mdb = totalLen / (2*bps); % a branch point is located on two different filaments, so technically they get counted twice. We divide by 2 to fix this
end