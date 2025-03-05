%% a program that launches Dragon and plots some results
close all
clear all
files = {'./Examples_ActinHypocotlys/kinza/bending3-new_curv/netProps.mat'};
% files = {'./Examples_ActinHypocotlys/kinza/bending3-new_curv/netProps.mat', ...
%     './Examples_ActinHypocotlys/kinza/classic3-new_curv/netProps.mat'};
% files = {'./Examples_ActinHypocotlys/kinza/bending3/netProps.mat', ...
%     './Examples_ActinHypocotlys/kinza/bending3-higher_lambda/netProps.mat', ...
%     './Examples_ActinHypocotlys/kinza/bending3-lam12/netProps.mat'};
% files = {'./Examples_ActinHypocotlys/kinza/test_classic/results_initial_param.mat', ...
%     './Examples_ActinHypocotlys/kinza/bending/results_bending.mat'};
%mylegend = {'Classic','Bending'} ;
%mylegend = {'Classic'} ;
mylegend = {'bending', 'classic'} ;
thr_length = 80 ; % minimal length to remove small hairs from the data
nbins = 20 ; % number of bins in histogram
markers = ['o', 's', 'd', '^', 'v', '>', '<', 'p', 'h'];

for i=1:length(files)
    
    % load data
    data = load(files{i}) ;
    results = data.net_props ;
    
    % extract data 
    length2D = results.filLenXY ; 
    curv_signed = results.curvatureSigned ;
    deviation = results.deviation ;
    curvature = results.newcurv ;
    
    % prepare filter
    filter = find(length2D>thr_length) ; % data to be kept - rest is removed
    
    % plot length
    figure(1);
    hold on
    histogram(length2D(filter) ,'BinWidth', 20, 'Normalization', 'probability'); % Adjust BinWidth as needed
    title('Histogram of Lengths');
    xlabel('Length in pix');
    ylabel('proba');
    legend(mylegend);
    hold off
    
    % plot curvature
    figure(2);
    hold on
    histogram(curv_signed(filter) ,'BinWidth', 0.005, 'Normalization', 'probability'); % Adjust BinWidth as needed
    title('Histogram of averaged signed curvature');
    xlabel('Averaged signed curvature in 1/pix');
    ylabel('proba');
    legend(mylegend);
    hold off
    
    % plot deviation
    figure(3);
    hold on
    histogram(deviation(filter), 'BinWidth', 1, 'Normalization', 'probability'); % Adjust BinWidth as needed
    title('Histogram of averaged deviation');
    xlabel('Averaged deviation in pix');
    ylabel('proba');
    xlim([0 40])
    legend(mylegend);
    hold off
    
    % plot curvature
    figure(4);
    hold on
    histogram(curvature(filter), 'BinWidth', 0.01, 'Normalization', 'probability'); % Adjust BinWidth as needed
    title('Histogram of averaged curvature');
    xlabel('Averaged curvature in 1/pix');
    ylabel('proba');
    xlim([0 1])
    legend(mylegend);
    hold off
    
    % create skeleton visualisation
    figure(4+i);
    skel = results.skelLabel ;
    num_labels = double(max(skel(:)));
    cmap = zeros(num_labels + 1, 3); % +1 for the background (0)
    cmap(1, :) = [0 0 0];  % Black
    cmap(2:end, :) = colormap(hsv(num_labels));
    hold on; % Important: Keep the plot open
    for k = 1:num_labels
        % filter too small length to be consistent with the rest of the
        % procedure
        if (curvature(k)>10)
            % Extract coordinates for the current label:
            [row, col] = find(skel == k);

            % Get the color for the current label from the colormap:
            color = cmap(k + 1, :); % +1 because the first color is for background (0)
            
            % Cycle through markers
            %marker = markers(mod(k-1, length(markers)) + 1); 
            
            % Plot the skeleton pixels with the corresponding color and line width:
            plot(col, row, '.', 'Color', color, 'MarkerSize', 3); % Adjust LineWidth as needed
    %         scatter(col, row, 'Marker', '.',...
    %                 'MarkerFaceColor', color, ...
    %                 'MarkerEdgeColor', color, 'HandleVisibility','off');
        end
    end
    title(mylegend{i})
    hold off
    axis equal;
    axis ij;
    
%     figure(4+length(files)+1); % attempt of violin plot
%     [f_normal, xi_normal] = ksdensity(deviation);
%     hold on
% 
%     % Plot for normal treatment
%     fill([xi_normal; flipud(xi_normal)], [f_normal; flipud(-f_normal)], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end