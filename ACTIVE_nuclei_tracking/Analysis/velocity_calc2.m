function [xyzs_id] = velocity_calc2(xyzs_id, xyzs_id_columns, accuracy)
% Function that uses finite difference theorom to calculate velocities
%
%  2/21/2013
%  R. Baker, M. Brasch
%
%  INPUTS:
%  xyzs_id: matrix of tracked cell information for all frames after
%           post-processing
%  xyzs_id_columns: column number containing cell IDs
%  accuracy: number of points to average over for finite difference theorom
%
%  OUTPUTS:
%  xyzs_id: matrix of tracked cell information after adding x and y
%           velocities as two additional columns
%
vx_column = size(xyzs_id,2) + 1;
vy_column = vx_column + 1;

% Sort data by cell ID and frame
xyzs_id = sortrows(xyzs_id, [xyzs_id_columns xyzs_id_columns-1]);

if nargin < 3;
   accuracy = 6; % Must be 2, 4, 6, or 8
end

% Determine the finite difference coefficients
switch accuracy
    case 2
        coefficient_vec = [-1/2 0 1/2];
    case 4
        coefficient_vec = [1/12 -2/3 0 2/3 -1/12];
    case 6
        coefficient_vec = [-1/60 3/20 -3/4 0 3/4 -3/20 1/60];
    case 8
        coefficient_vec = [1/280 -4/105 1/5 -4/5 0 4/5 -1/5 4/105 -1/280];
end

zero_index = ceil(length(coefficient_vec)/2); % Find the center of the tridiagonal
xyzs_id_index = (1:size(xyzs_id,1))';
num_cells = max(xyzs_id(:,xyzs_id_columns));

for i=1:num_cells;
    % Extract information for the given cell
    cell_bool = xyzs_id(:,xyzs_id_columns)==i;
    xpos_temp = xyzs_id(cell_bool,1);
    ypos_temp = xyzs_id(cell_bool,2);
    frame_temp = xyzs_id(cell_bool,xyzs_id_columns-1);
    index_temp = xyzs_id_index(cell_bool);
    velocity_keep_vec = ones(1,length(frame_temp));
    
    % Check for missing frames, which will be removed from analysis
    for t = 2:length(frame_temp);        
        if frame_temp(t) ~= frame_temp(t-1) + 1;
            if t <= (accuracy/2) + 1
                velocity_keep_vec(1:t+(accuracy/2)-1) = 0;
            elseif t >= length(frame_temp)-(accuracy/2)+1;
                velocity_keep_vec(t-(accuracy/2):length(frame_temp)) = 0;
            else
                velocity_keep_vec(t-(accuracy/2):t+(accuracy/2)-1) = 0;
            end
        end
    end
    
    % Positions near the beginning and end will not be counted
    if length(frame_temp) <= accuracy/2
        velocity_keep_vec = velocity_keep_vec*0;
    else
        velocity_keep_vec(1:accuracy/2)=0;
        velocity_keep_vec((length(frame_temp)-accuracy/2+1):length(frame_temp)) = 0;
    end
    
     if length(frame_temp) >= length(coefficient_vec);
        n_points = length(frame_temp);
        diag_matrix = zeros(n_points,n_points); % Initialize the size of the coefficient matrix
        
        for k = 1:length(coefficient_vec) % loop that creates the coefficient matrix
            off_diag_size = ones(n_points - (abs(zero_index-k)),1); % Determine the size of one of the off diagonals
            diag_matrix = diag_matrix + diag(coefficient_vec(k)*off_diag_size,-(zero_index-k)); % Update the diagonal matrix
        end
        % Calculate the velocities in x and y using finite difference
        % theorom
        vx_vector = diag_matrix*xpos_temp;
        vy_vector = diag_matrix*ypos_temp;
        
        % Get rid of velocities affected by missing frames
        vx_vector = vx_vector.*velocity_keep_vec';
        vy_vector = vy_vector.*velocity_keep_vec';
        
        % Insert velocity data into xyzs_id matrix
        for z=1:length(frame_temp)
            xyzs_id(index_temp(z),vx_column) = vx_vector(z);
            xyzs_id(index_temp(z),vy_column) = vy_vector(z);
        end
    end
end


