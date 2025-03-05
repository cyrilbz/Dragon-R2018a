function b = getBinary(BW, msk, prc_thresh)
% GETBINARY  Produces binary image from greyscale image of filamentous content
%
%   @input: BW - greyscale image of filamentous content
%           msk - binary mask of the region of interest
%           prc_thresh - percentile value for thresholding
%
%   @output: b - binary image

    im_fib = fibermetric(BW);
    msk = msk > 0; %ensure binary
    
    %%%% version withou stat and machnie learning toolbox
    im_fib_mask = im_fib(msk) ; % mask image
    sorted = sort(im_fib_mask(:)); % sort values
    n_pxl = length(sorted); % number of pixels
    thresh_idx = max(1, round(prc_thresh/100*n_pxl)) ; % find the indice for the given percentile value
    thresh_value = sorted(thresh_idx) ; % get threshold value
    b = im_fib > thresh_value ;
    %%%% previous version that requires Stat and Machine learning toolbox
    %%% b = im_fib > prctile(im_fib(msk),prc_thresh,'all'); 
    b = bwareaopen(b,50); %remove very small fragments
    b(~msk) = 0;
end