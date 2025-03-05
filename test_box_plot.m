% Example data
normal_length = [1.2, 1.5, 1.3, 1.4];
bending_length = [1.6, 1.7, 1.5, 1.8];

% Compute quartiles for each dataset
q_normal = quartiles(normal_length);
q_bending = quartiles(bending_length);

% Define positions for plotting
positions = [1, 2];

% Create box plot
figure;
hold on;

% Plot boxes
box_width = 0.2;
fill([positions(1)-box_width, positions(1)+box_width, positions(1)+box_width, positions(1)-box_width], ...
     [q_normal(1), q_normal(1), q_normal(3), q_normal(3)], 'b', 'FaceAlpha', 0.5);
fill([positions(2)-box_width, positions(2)+box_width, positions(2)+box_width, positions(2)-box_width], ...
     [q_bending(1), q_bending(1), q_bending(3), q_bending(3)], 'r', 'FaceAlpha', 0.5);

% Plot medians
plot([positions(1)-box_width, positions(1)+box_width], [q_normal(2), q_normal(2)], 'b-');
plot([positions(2)-box_width, positions(2)+box_width], [q_bending(2), q_bending(2)], 'r-');

% Plot whiskers
plot([positions(1), positions(1)], [min(normal_length), q_normal(1)], 'b-');
plot([positions(1), positions(1)], [q_normal(3), max(normal_length)], 'b-');
plot([positions(2), positions(2)], [min(bending_length), q_bending(1)], 'r-');
plot([positions(2), positions(2)], [q_bending(3), max(bending_length)], 'r-');

% Add labels and title
set(gca, 'XTick', positions);
set(gca, 'XTickLabel', {'Normal', 'Bending'});
ylabel('Length');
title('Box Plot of Length by Treatment');

hold off;
% Function to compute quartiles
function q = quartiles(data)
    sorted_data = sort(data);
    n_val = length(sorted_data); % number of values
    max_idx = max(1, round(0.95 * n_val)); % index for 95th percentile
    min_idx = max(1, round(0.05 * n_val)); % index for 5th percentile
    q1 = sorted_data(max_idx); % 95th percentile value
    q3 = sorted_data(min_idx); % 5th percentile value
    q1 = median(sorted_data(1:floor(n_val/2)));
    q3 = median(sorted_data(ceil(n_val/2):end));
    q = [q1, median(sorted_data), q3];
end