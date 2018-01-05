function [lambda1, lambda2, radius_gyration, track_length, gyration_angle, asphericity] = gyration(xyzs_id,xyzs_id_columns)
% Function that determines the parameters for the gyration tensor of cell
% tracks
%
%  6/20/2013
%  R. Baker, M. Brasch
%
%  INPUTS:
%  xyzs_id: matrix of tracked cell information for all frames after
%           post-processing
%  xyzs_id_columns: column number containing cell IDs
%
%  OUTPUTS:
%  lambda1, lambda2: vector containing amsllest and largest principle 
%           moments of gyration tensor, respectively, for each cell
%  radius_gyration: vector of the radius of gyration for each cell
%  track_length: vector of track lengths for each cell
%  gyration_angle: vector of principle axis (largest eigenvector) for each
%           cell
%  asphericity: vector of cell track asphericities as determined from 

% Sort by cell ID and frame
xyzs_id = sortrows(xyzs_id,[xyzs_id_columns xyzs_id_columns-1]);
% Generate list of unique cell IDs
cell_IDs = unique(xyzs_id(:,xyzs_id_columns));

% Pre-allocate gyration parameters
gyration_tensor = cell(size(cell_IDs,1),1);
eigen_values = cell(size(cell_IDs,1),1);
eigen_vectors = cell(size(cell_IDs,1),1);
track_length = zeros(size(cell_IDs,1),1);
lambda1 = zeros(size(cell_IDs,1),1);
lambda2 = zeros(size(cell_IDs,1),1);
radius_gyration = zeros(size(cell_IDs,1),1);
gyration_angle = zeros(size(cell_IDs,1),1);
asphericity = zeros(size(cell_IDs,1),1);

% Loop that calculates gyration tensory, principle moments, and shape
% descriptors of cell tracks
for a=1:size(cell_IDs,1);
    bool = xyzs_id(:,xyzs_id_columns) == cell_IDs(a);
    cell_mat = xyzs_id(bool,:);
    
    if size(cell_mat,1) < 3
        continue
    end
    
    Sxy = 0;
    Sxx = 0;
    Syy = 0;
    % Calculate each component of gyration tensor
    for i=1:size(cell_mat,1);
        for j=1:size(cell_mat,1);
            Sxy = Sxy + (cell_mat(i,1)-cell_mat(j,1))*(cell_mat(i,2)-cell_mat(j,2));
            Sxx = Sxx + (cell_mat(i,1)-cell_mat(j,1))*(cell_mat(i,1)-cell_mat(j,1));
            Syy = Syy + (cell_mat(i,2)-cell_mat(j,2))*(cell_mat(i,2)-cell_mat(j,2));
        end
    end
    Sxy = Sxy/(2*(size(cell_mat,1))^2);
    Sxx = Sxx/(2*(size(cell_mat,1))^2);
    Syy = Syy/(2*(size(cell_mat,1))^2);
  
    gyration_tensor{a} = [Sxx Sxy; Sxy Syy];
  
    % Calculate principle moments
    [eigen_vectors{a} eigen_values{a}] = eig(gyration_tensor{a});
    
    % Calculate angle of rotation for largest principle moment
    gyration_angle(a) = atan2(eigen_vectors{a}(2,2),eigen_vectors{a}(1,2));

    % Determine gyration tensor parameters and shape descriptors
    lambda1(a) = eigen_values{a}(1,1);
    lambda2(a) = eigen_values{a}(2,2);
    radius_gyration(a) = lambda1(a) + lambda2(a);
    track_length(a) = cell_mat(j,xyzs_id_columns-1)-cell_mat(1,xyzs_id_columns-1);
    asphericity(a) = (lambda2(a)-lambda1(a))./(lambda2(a)+lambda1(a));
end

% Remove data for tracks with less than 2 points
lambda1(radius_gyration==0)=[];
lambda2(radius_gyration==0)=[];
track_length(radius_gyration==0)=[];
asphericity(radius_gyration==0)=[];
radius_gyration(radius_gyration==0)=[];
end

