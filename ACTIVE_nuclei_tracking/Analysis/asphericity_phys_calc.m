function [asphericity_phys, asphericity_phys_avg] = asphericity_phys_calc(lambda1, lambda2, track_length)
% Function that calculates track asphericities
%
%  2/8/2012
%  M. Brasch, R. Baker, L. Manning
%
%  INPUTS:
%  lambda1, lambda2: cell_arrays containing amsllest and largest principle 
%           moments of gyration tensor, respectively, for each replicate
%  track_length: vector of track lengths for each cell 
%
%  OUTPUTS:
%  asphericity_phys: cell array contain track asphericity for every cell 
%  asphericity_phys_avg: 
%
% LOCAL PARAMETERS:
track_length_cutoff = 0; % Minimum track size to be counted in analysis

% Preallocation
asphericity_phys = cell(size(lambda1,1),1); % Will contain every track's asphericity data based on the physics equation
asphericity_phys_avg = zeros(size(lambda1,1),2); % Will contain avg and #of tracks 

% 
for i=1:size(lambda1,1)
    % Calculate every track's asphericity
    bool = track_length{i}(:,1) >= track_length_cutoff;
    asphericity_phys{i} = (lambda2{i}(bool)-lambda1{i}(bool))./(lambda2{i}(bool)+lambda1{i}(bool));    

    % Calculate the average and #of tracks for each replicate
    asphericity_phys_avg(i,1) = mean(asphericity_phys{i}(:,1));
    asphericity_phys_avg(i,2) = size(lambda1{i},1);
end

% Group the replicates together to calculate and average, #of tracks, and
% std (of the replicate means)
n_groups = 9; % Number of total groups
n_groups1 = 3; % Number of substrates
n_groups2 = 3; % Number of densities
n_samples = 4; % Number of replicates
index = 0;
group_asphericity = zeros(n_groups,3);
for r=1:n_groups
    for s=1:n_samples
       index = index+1;
       group_asphericity(r,1) = (group_asphericity(r,1)*group_asphericity(r,2)+asphericity_phys_avg(index,1)*asphericity_phys_avg(index,2))/(asphericity_phys_avg(index,2)+group_asphericity(r,2));
       group_asphericity(r,2) = group_asphericity(r,2)+asphericity_phys_avg(index,2);
    end
    % Replicate std (std of replicate means)
    group_asphericity(r,3) = std(asphericity_phys_avg(index-n_samples+1:index));
end

% Plot the group average asphericity and standard deviation
index = 1;
color{1} = 'b'; color{2} = 'r'; color{3} = 'g';
for i=1:n_groups1;
    for j=1:n_groups2;
        plot_mat1(i,j) = group_asphericity(index,1);
        error_mat1(i,j) = group_asphericity(index,3);
        index = index + 1;
    end
    
    figure(800)
    errorbar([3 1 2],plot_mat1(i,:),error_mat1(i,:),[color{i},'o']);
    hold on
end

% Set Figure Parameters
figure(800)
labels = {'5k cells/cm^2' '10k cells/cm^2' '20k cells/cm^2'};
set(gca, 'XTick', 1:3, 'XTickLabel', labels);
ylabel('Normalized Asphericity (pixels^2/frame)')
legend('Wrinkled', 'Non-wrinkled', 'TCPS','Location', 'EastOutside')
title('Summarized Asphericity Data')

end
