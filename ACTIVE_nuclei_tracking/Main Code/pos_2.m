function [ lookup_table, pos_matrix ] = pos_2( event_info, xyzs_id, xyzs_id_columns, frame_avg )
%   The function performs a position cost analysis for each non-duplicated 
%   entry in a specific event's information. This information is used to
%   correct for mislabelled cells during tracking
% 
%  2/8/2012
%  M. Brasch, R. Baker, L. Manning
%   
% Conditions and terms of use:
% The software packages provided here are M-files executable in MATLAB, a 
% proprietary numerical computing enviornment developed by MathWorks.
% You are free to use this software for research purposes, but you should 
% not redistribute it without the consent of the authors. In addition, end 
% users are expected to include adequate citations and acknowledgments 
% whenever results or derivatives that are based on the software are presented or published.
%
% Citation to ACTIVE should include the following:
% Baker RM, Brasch ME, Manning ML, and Henderson JH. Automated, 
%        contour-based tracking and analysis of cell behavior over long 
%        timescales in environments of varying complexity and cell density.
%        Journal information to be updated when available.
%
% Citations to work foundational to ACTIVE are suggested to include the following, at a minimum:
%
% Idema T. A new way of tracking motion, shape, and divisions. European 
%        Biophysics Journal. 2013:1-8.
% Crocker JC, Grier DG. Methods of digital video microscopy for colloidal 
%        studies. Journal of Colloid and Interface Science. 1996;179(1):298-310.
% Gao Y, Kilfoil ML. Accurate detection and complete tracking of large 
%        populations of features in three dimensions. Optics Express. 
%        2009;17(6):4685-704.
%
%  INPUTS:
%  event_info: matrix of particle information (from collision_corrector.m)
%              for a given event
%  xyzs_id: matrix of particle information post-tracking (output from
%           trackmem_new.m)
%  frame_avg: contains average values for cell major axis, minor axis, 
%             aspect ratio, angle of orientation, area, intensity,
%             integrated intensity
%  OUTPUTS:
%  lookup_table: 4x4 matrix of information from the pos_matrix of both
%                cells in an event in the first and last frames of the event; used to
%                determine if cell IDs should be switched
%  pos_matrix: matrix containing position analysis for all frames of an
%              event (an output for debugging purposes only)
%
%  LOCAL PARAMETERS:

%Isolate event cell information
bool_pos = (event_info(:,6) == 6) | (event_info(:,6) == 7);
num_dup = sum(bool_pos);
pos_size = 2*(size(event_info,1) - num_dup);
pos_matrix = zeros(pos_size,10);
switch_vec = zeros(pos_size/2,1);
pos_count = 1;              %Counter for pos_matrix
cell1 = event_info(1,2);
cell2 = event_info(1,3);
cell_index_vector = zeros(pos_size,1);

max_consec_dup = 3;
dup_count = 0;


%Organize position matrix for all cell ID entries in event info
for i = 1:size(event_info,1)
    
    if event_info(i,6) == 6 || event_info(i,6) == 7
        dup_count = dup_count + 1;
        continue
    end
    
        cell1_index = find((xyzs_id(:,xyzs_id_columns) == event_info(i,2)) & (xyzs_id(:,xyzs_id_columns-1) == event_info(i,4)));
        cell2_index = find((xyzs_id(:,xyzs_id_columns) == event_info(i,3)) & (xyzs_id(:,xyzs_id_columns-1) == event_info(i,4)));
        pos_matrix(pos_count,1) = xyzs_id(cell1_index, xyzs_id_columns - 1); % frame 1
        pos_matrix(pos_count,2) = xyzs_id(cell1_index,1); % x coord cell 1
        pos_matrix(pos_count,3) = xyzs_id(cell1_index,2); % y coord cell 1 
        pos_matrix(pos_count,4) = xyzs_id(cell1_index,xyzs_id_columns); % cell 1 ID
        pos_matrix(pos_count,9) = event_info(i,1);
        pos_matrix(pos_count,11) = cell1_index;
        pos_matrix(pos_count+1,1) = xyzs_id(cell2_index,xyzs_id_columns - 1); % frame 1
        pos_matrix(pos_count+1,2) = xyzs_id(cell2_index,1); % x coord cell 2
        pos_matrix(pos_count+1,3) = xyzs_id(cell2_index,2); % y coord cell 2
        pos_matrix(pos_count+1,4) = xyzs_id(cell2_index,xyzs_id_columns); % cell ID 2
        pos_matrix(pos_count+1,9) = event_info(i,9);
        pos_matrix(pos_count+1,11) = cell2_index;
        
        cell_index_vector(pos_count) = cell1_index;
        cell_index_vector(pos_count+1) = cell2_index;
        
        if dup_count >= max_consec_dup
            pos_matrix(pos_count,10) = 1;
            pos_matrix(pos_count+1,10) = 1;
        end    
        
        dup_count = 0;
        pos_count = pos_count+2;    
end

%Initialize first frame information for both cells
pos_matrix(1,5) = cell1;
pos_matrix(2,5) = cell2;

%The same counter is used for inserting position switch information
pos_count = 3;

pos_matrix(1,8) = cell1;
pos_matrix(2,8) = cell2;


while pos_count < size(pos_matrix,1)
    % if number of consecutive duplicates is more than max allowed, use
    % fingerprint method for comparison, otherwise use position method
    if pos_matrix(pos_count,10) == 1
        [cost13,~] = cost(xyzs_id(pos_matrix(pos_count-2,11),:),xyzs_id(pos_matrix(pos_count,11),:), frame_avg, xyzs_id_columns);
        [cost24,~] = cost(xyzs_id(pos_matrix(pos_count-1,11),:),xyzs_id(pos_matrix(pos_count+1,11),:), frame_avg, xyzs_id_columns);
        cost1 = cost13 + cost24;
        
        [cost14,~] = cost(xyzs_id(pos_matrix(pos_count-2,11),:),xyzs_id(pos_matrix(pos_count+1,11),:), frame_avg, xyzs_id_columns);
        [cost23,~] = cost(xyzs_id(pos_matrix(pos_count-1,11),:),xyzs_id(pos_matrix(pos_count,11),:), frame_avg, xyzs_id_columns);
        cost2 = cost14 + cost23;
    else
        % Calculate the distance between centroids
        cost13 = ((pos_matrix(pos_count,2) - pos_matrix(pos_count-2,2))^2 + (pos_matrix(pos_count,3) - pos_matrix(pos_count-2,3))^2)^(1/2); 
        cost24 = ((pos_matrix(pos_count+1,2) - pos_matrix(pos_count-1,2))^2 + (pos_matrix(pos_count+1,3) - pos_matrix(pos_count-1,3))^2)^(1/2); 
        cost1 = cost13 + cost24;

        cost14 = ((pos_matrix(pos_count+1,2) - pos_matrix(pos_count-2,2))^2 + (pos_matrix(pos_count+1,3) - pos_matrix(pos_count-2,3))^2)^(1/2); 
        cost23 = ((pos_matrix(pos_count,2) - pos_matrix(pos_count-1,2))^2 + (pos_matrix(pos_count,3) - pos_matrix(pos_count-1,3))^2)^(1/2); 
        cost2 = cost14+cost23;
    end

    pos_matrix(pos_count,6) = cost1;
    pos_matrix(pos_count+1,6) = cost2;

    if cost1 < cost2
        pos_matrix(pos_count,5) = pos_matrix(pos_count-2,5);
        pos_matrix(pos_count+1,5) = pos_matrix(pos_count-1,5);
    else 
        pos_matrix(pos_count,5) = pos_matrix(pos_count-1,5);
        pos_matrix(pos_count+1,5) = pos_matrix(pos_count-2,5);
    end
    
    %Advance position counter
    pos_count = pos_count+2;
end

%Create lookup_table (will be reorganized and output into mastertable)
lookup_table(1,1) = pos_matrix(1,1);
lookup_table(1,2) = pos_matrix(1,4);
lookup_table(1,3) = pos_matrix(1,5);
lookup_table(1,4) = pos_matrix(1,9);
lookup_table(2,1) = pos_matrix(2,1);
lookup_table(2,2) = pos_matrix(2,4);
lookup_table(2,3) = pos_matrix(2,5);
lookup_table(2,4) = pos_matrix(2,9);
lookup_table(3,1) = pos_matrix(pos_count-2,1);
lookup_table(3,2) = pos_matrix(pos_count-2,4);
lookup_table(3,3) = pos_matrix(pos_count-2,5);
lookup_table(3,4) = pos_matrix(pos_count-2,9);
lookup_table(4,1) = pos_matrix(pos_count-1,1);
lookup_table(4,2) = pos_matrix(pos_count-1,4);
lookup_table(4,3) = pos_matrix(pos_count-1,5);
lookup_table(4,4) = pos_matrix(pos_count-1,9);

end