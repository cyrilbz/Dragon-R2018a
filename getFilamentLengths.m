function len = getFilamentLengths(lbl)
% GETFILAMENTLENGTHS  determines the length of filaments from LABELMATRIX
%
%   @input: lbl - a labelled matrix
%
%   @output: len - array containing filament lengths; index corresponds
%   to the label number from lbl
%
%   Any lengths of NaN returned would suggest that index doesn't exist.
%   This is done so that mead and std calculations remain accurate with the
%   correct flag options
%
%   If one labelled filament (8-connected) becomes n connected filaments
%   when 4-connected, then there are n-1 diagonal joins, each adding 
%   (sqrt(2)-1) to the length.

    no_of_filaments = max(lbl(:));
    len = NaN(no_of_filaments,1);
    extra_diag_len = sqrt(2)-1; %moving diagonally is longer by this much
    for i = 1:no_of_filaments
        filament = (lbl == i); %grab a single filament in the entire image

        % to check results using euclidean norm
%         endpoints = (bwmorph(filament, 'endpoints'));
%         [row, col] = find(endpoints);
%         endpt=[row(1),col(1)];
%         startpt=[row(2),col(2)];
%         disp(norm(endpt-startpt))
        
        % number of pixels in the filament
        px_sum = sum(filament(:),'omitnan');

        % Identify 4-connected pixels
        sk_structure_4_conn = bwareaopen(filament, 2, 4);
        sum_sk_structure_4_conn = sum(sk_structure_4_conn(:));
        
        % Compute number of 8-connected pixels
        sum_sk_structure_8_conn = px_sum - sum_sk_structure_4_conn;

        % Compute length multiplying diagonal pixels with sqrt(2)
        len(i) = (sum_sk_structure_4_conn * 1) +...
            (sum_sk_structure_8_conn * 1 * sqrt(2));

          % old computation that generated an error
%         filament_lbl = bwlabel(filament,4); %relabel with 4-connected rules
%         n = max(filament_lbl(:),'omitnan') - 1; %num of diagonal joins
%         len(i) = px_sum + (n*extra_diag_len);
    end
end