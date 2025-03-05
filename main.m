%% a program that launches Dragon and plots some results
close all
clear all

%results = Dragon("./Examples_ActinHypocotlys/Col0/001/") ;
%results = Dragon("./Examples_ActinHypocotlys/kinza/bending3-test_curv/") ;
results = Dragon("./Examples_ActinHypocotlys/kinza/bending-new_curv/") ;
results = Dragon("./Examples_ActinHypocotlys/kinza/bending2-new_curv/") ;
results = Dragon("./Examples_ActinHypocotlys/kinza/bending4-new_curv/") ;
%results2 = Dragon("./Examples_ActinHypocotlys/kinza/classic3-higher_lambda/") ;
% 
%% plots (use "run section" to avoid computing data agin!)
% close all
% origin = results.imRaw ;
% binary = results.imMask ;
% skel = results.skelLabel ;
% skeleton = results.imSkel ;
% figure(1);
% imshow(origin)
% 
% % figure(3);
% % imshow(skeleton)
% 
% % show skeleton
% figure(2);
% num_labels = double(max(skel(:)));
% cmap = zeros(num_labels + 1, 3); % +1 for the background (0)
% cmap(1, :) = [0 0 0];  % Black
% cmap(2:end, :) = colormap(jet(num_labels));
% % imshow(skel,[0 num_labels]); % Important: Set display range!
% % colormap(cmap);
% hold on; % Important: Keep the plot open
% for i = 1:num_labels
%     % Extract coordinates for the current label:
%     [row, col] = find(skel == i);
% 
%     % Get the color for the current label from the colormap:
%     color = cmap(i + 1, :); % +1 because the first color is for background (0)
% 
%     % Plot the skeleton pixels with the corresponding color and line width:
% %     plot(col, row, '.', 'Color', color, 'LineWidth', 2); % Adjust LineWidth as needed
%     scatter(col, row, 'Marker', '.',...
%             'MarkerFaceColor', color, ...
%             'MarkerEdgeColor', color, 'HandleVisibility','off');
% end
% hold off
% axis equal;
% axis ij;
% 
% % plot length histogram
% 
% l_max = 80 ; % max length to take into account statistics
% length2D = results.filLenXY ;
% filter = find(length2D>l_max) ; % data to be kept - rest is removed
% 
% figure(4);
% histogram(length2D(filter) , 20, 'Normalization', 'probability'); % Adjust BinWidth as needed
% title('Histogram of Lengths');
% xlabel('Length in pix');
% ylabel('Frequency');
% 
% % plot NEW curvature histogram
% curv = results.newcurv ;
% figure(5);
% histogram(curv(filter) ,40, 'Normalization', 'probability'); % Adjust BinWidth as needed
% title('Histogram of curvatures');
% xlabel('curv');
% ylabel('Frequency');
% 
% % plot Signed curvature histogram
% curv_signed = results.curvatureSigned ;
% figure(6);
% histogram(curv_signed(filter) ,20, 'Normalization', 'probability'); % Adjust BinWidth as needed
% title('Histogram of signed curvatures');
% xlabel('curv');
% ylabel('Frequency');
% 
% % plot deviation histogram
% deviation = results.deviation ;
% figure(7);
% histogram(deviation(filter) ,20, 'Normalization', 'probability'); % Adjust BinWidth as needed
% title('Histogram of deviation');
% xlabel('deviation');
% ylabel('Frequency');