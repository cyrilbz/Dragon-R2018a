clear all
close all
set(0, 'DefaultFigureColor', 'w');
%% a prog to aggregate results and compute statistics

% list data
type1 = {'./Examples_ActinHypocotlys/kinza/bending-new_curv/netProps.mat', ...
    './Examples_ActinHypocotlys/kinza/bending2-new_curv/netProps.mat', ...
    './Examples_ActinHypocotlys/kinza/bending3-new_curv/netProps.mat', ...
    './Examples_ActinHypocotlys/kinza/bending4-new_curv/netProps.mat'};

type2 = {'./Examples_ActinHypocotlys/kinza/classic-new_curv/netProps.mat', ...
    './Examples_ActinHypocotlys/kinza/classic2-new_curv/netProps.mat', ...
    './Examples_ActinHypocotlys/kinza/classic3-new_curv/netProps.mat', ...
    './Examples_ActinHypocotlys/kinza/classic4-new_curv/netProps.mat'};

% type1 = {'./Examples_ActinHypocotlys/kinza/bending3-new_curv/netProps.mat'};
% 
% type2 = {'./Examples_ActinHypocotlys/kinza/classic3-new_curv/netProps.mat'};

mylegend = {'Bending 35 min(n=4)','Straight (n=4)'} ;

thr_length = 80 ; % minimal length to remove small hairs from the data
nbins = 20 ; % number of bins in histogram
markers = ['o', 's', 'd', '^', 'v', '>', '<', 'p', 'h'];

%% aggregate type 1
length_type1 = [] ;
curvature_type1 = [] ;
deviation_type1 = [] ;
density_type1 = [] ;
for i=1:length(type1)
    % load data
    data = load(type1{i}) ;
    results = data.net_props ;
    
    % extract and aggregate
    length_type1 = [length_type1 ;  results.filLenXY] ; 
    curvature_type1 = [curvature_type1 , results.newcurv];
    deviation_type1 = [deviation_type1 , results.deviation];
    density_type1 = [density_type1, results.skelDensity] ;
end 

%% aggregate type 2
length_type2 = [] ;
curvature_type2 = [] ;
deviation_type2 = [] ;
density_type2 = [] ;
for i=1:length(type2)
    % load data
    data = load(type2{i}) ;
    results = data.net_props ;
    
    % extract and aggregate
    length_type2 = [length_type2 ;  results.filLenXY] ; 
    curvature_type2 = [curvature_type2 , results.newcurv];
    deviation_type2 = [deviation_type2 , results.deviation];
    density_type2 = [density_type2, results.skelDensity] ;
end 
fprintf('Mean density 1 : %.4f\n', mean(density_type1));
fprintf('Mean density 2 : %.4f\n', mean(density_type2));

%% plots

% prepare filters
filter1 = find(length_type1>thr_length) ; % data to be kept - rest is removed
filter2 = find(length_type2>thr_length) ; % data to be kept - rest is removed

% plot length histograms
figure(1);
hold on
histogram(length_type1(filter1) ,'BinWidth', 60, 'Normalization', 'probability'); 
histogram(length_type2(filter2) ,'BinWidth', 60, 'Normalization', 'probability'); 
% histogram(length_type1(filter1) , 25, 'Normalization', 'probability'); 
% histogram(length_type2(filter2) , 25, 'Normalization', 'probability'); 
title('Lengths');
xlabel('Length in pix');
ylabel('proba');
legend(mylegend);
hold off

% plot curvature
figure(2);
hold on
histogram(curvature_type1(filter1) ,'BinWidth', 0.01, 'Normalization', 'probability'); 
histogram(curvature_type2(filter2) ,'BinWidth', 0.01, 'Normalization', 'probability'); 
% histogram(curvature_type1(filter1) ,25, 'Normalization', 'probability'); 
% histogram(curvature_type2(filter2) ,25, 'Normalization', 'probability'); 
title('Averaged curvature');
xlabel('Averaged curvature in 1/pix');
ylabel('proba');
xlim([0 2])
legend(mylegend);
hold off

% plot deviation
figure(3);
hold on
histogram(deviation_type1(filter1) ,'BinWidth', 1, 'Normalization', 'probability'); 
histogram(deviation_type2(filter2) ,'BinWidth', 1, 'Normalization', 'probability'); 
% histogram(curvature_type1(filter1) ,25, 'Normalization', 'probability'); 
% histogram(curvature_type2(filter2) ,25, 'Normalization', 'probability'); 
title('Averaged deviation');
xlabel('Averaged deviation in pix');
ylabel('proba');
xlim([0 40])
legend(mylegend);
hold off


%% box plots

boxplot(deviation_type1(filter1),deviation_type2(filter2),'deviation',mylegend)
boxplot(curvature_type1(filter1),curvature_type2(filter2),'curvature',mylegend)
boxplot(length_type1(filter1),length_type2(filter2),'length',mylegend)


function boxplot(data1,data2,name,mylegend)
% data
stat1 = calculateStatistics(data1);
stat2 = calculateStatistics(data2);
% normal_length = data1 ;
% bending_length = data2 ;

% Compute quartiles for each dataset
q_1 = quartiles(data1);
q_2 = quartiles(data2);

% Define positions for plotting
positions = [1, 2];

% Create box plot
figure;
hold on;

% Plot boxes
box_width = 0.2;
fill([positions(1)-box_width, positions(1)+box_width, positions(1)+box_width, positions(1)-box_width], ...
     [q_1(1), q_1(1), q_1(3), q_1(3)], 'b', 'FaceAlpha', 0.5);
fill([positions(2)-box_width, positions(2)+box_width, positions(2)+box_width, positions(2)-box_width], ...
     [q_2(1), q_2(1), q_2(3), q_2(3)], 'r', 'FaceAlpha', 0.5);

% Plot medians
plot([positions(1)-box_width, positions(1)+box_width], [q_1(2), q_1(2)], 'b-');
plot([positions(2)-box_width, positions(2)+box_width], [q_2(2), q_2(2)], 'r-');

% Plot means
plot([positions(1)-box_width, positions(1)+box_width], [stat1.mean, stat1.mean], 'b--');
plot([positions(2)-box_width, positions(2)+box_width], [stat2.mean, stat2.mean], 'r--');

% Plot whiskers
plot([positions(1), positions(1)], [stat1.min_value, q_1(1)], 'b-');
plot([positions(1), positions(1)], [q_1(3), stat2.max_value], 'b-');
plot([positions(2), positions(2)], [stat2.min_value, q_2(1)], 'r-');
plot([positions(2), positions(2)], [q_2(3), stat2.max_value], 'r-');

% Add labels and title
set(gca, 'XTick', positions);
set(gca, 'XTickLabel', {mylegend{1}, mylegend{2}});
ylabel(name);
title(name);
hold off;
end

% Function to compute quartiles
function q = quartiles(data)
    sorted_data = sort(data);
    n_val = length(sorted_data); % number of values
    q1 = median(sorted_data(1:floor(n_val/2))); % 25% value
    q3 = median(sorted_data(ceil(n_val/2):end)); % 75% value
    q = [q1, median(sorted_data), q3];
end

% function to compute statistics
function stats = calculateStatistics(data)
    % Calculate mean
    stats.mean = mean(data);

    % Calculate median
    stats.median = median(data);

    % Calculate 95% interval
    sorted = sort(data);
    n_val = length(sorted); % number of values
    max_idx = max(1, round(0.95 * n_val)); % index for 95th percentile
    min_idx = max(1, round(0.05 * n_val)); % index for 5th percentile
    stats.max_value = sorted(max_idx); % 95th percentile value
    stats.min_value = sorted(min_idx); % 5th percentile value
    
%     % Display statistics
%     fprintf('Mean: %.2f\n', stats.mean);
%     fprintf('Median: %.2f\n', stats.median);
%     fprintf('95%% Interval: [%.2f, %.2f]\n', stats.min_value, stats.max_value);
end