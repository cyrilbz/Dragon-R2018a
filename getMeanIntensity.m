function in = getMeanIntensity(I, lbl, widths, filAng) %ans.filamentAng
% GETMEANINTENSITY  Gets the mean intensity of a filament using avg. width
%
%   @input: I      - image to grab intensities from (np.imRaw)
%           lbl    - labelled network of filaments, obtained from DRAGoN
%           widths - list of filament widths, from DRAGoN
%           filAng - structure(s) of filament angles, from DRAGoN
%                    (np.filamentAng as seen in getNetworkProperties.m)
%
%   @output: in - list of mean intensities of each labelled filament (index
%                 represents label number)
%
%   To be able to just look at the centreline intensity, remove the
%   imdilate function call and use: msk = (lbl==idx);
%
%   To call within getNetworkProperties.m, use the following:
%       np.meanIn = getMeanIntensity(np.imRaw, np.skelLabel, np.filWidth, np.filamentAng);

    total_fils = length(widths);
    in = zeros(total_fils,1);
    for idx = 1:total_fils
        ang = rad2deg(filAng(idx).angXY);
        SE = strel("line", widths(idx), ang);
        msk = imdilate(lbl==idx, SE);
        in(idx) = mean(I(msk), "all");
    end
end