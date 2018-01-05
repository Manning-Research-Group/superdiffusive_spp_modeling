function [xyzs_id, sib_matrix, event_array, xyzs_id_columns] = collision_corrector_rtc2(xyzs_id, sib_matrix, xyzs_id_columns,max_collision_time)
% Function to construct events between sibling cells, identify them as
% collisions or divisions, and then create an array containing all of that
% information to be further analyzed by collision_tags
% 
%  11/8/2012
%  R. Baker, M. Brasch, L. Manning
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
%  xyzs_id: matrix of tracked cell information for all frames (output from
%           trackmem_new)
%  sib_matrix: matrix of sibling information (output from sib_matrix_creation) 
%  xyzs_id_columns: the column containing cell IDs in the xyzs_id matrix
%  max_collision_time: max time a complete cell occlusion occurs
%
%  OUTPUTS:
%  sib_matrix: 4 new columns added to sib_matrix including 1)type of event,
%              2)index for a continuing event, 3)event number, and 4)index of the
%              second cell into xyzs_id
%  event_array: cell array containing the event information for each
%               event, including information for each cell for every frame of a given
%               event
%
%  LOCAL PARAMETERS:
if nargin < 4
    max_collision_time = 10; % Max time for a cell collision
end

% Add indexing to the sib_matrix 
% Column 6 is the type of event, column 7 is an index to a previous event
% (for continuations or exits) and column 8 is the event tag (grouping all
% events of the same overall event together)
sib_matrix(:,6:9) = zeros(size(sib_matrix(1)));
for p=1:size(sib_matrix,1)
    cell2_index = sib_matrix(p,3)==sib_matrix(:,2) & sib_matrix(:,4) == sib_matrix(p,4);
    if ~isempty(sib_matrix(cell2_index,1))
        sib_matrix(p,9) = sib_matrix(cell2_index,1);
    end
end
sib_bool = sib_matrix(:,2) < sib_matrix(:,3);
sib_matrix = sib_matrix(sib_bool,:);
sib_matrix = sortrows(sib_matrix,4);

% Initialize number of wrong divisions
wrong_division_count = 0;    
% loop to determine if it is a sibling, collision, or continuing collision
for i = 1:size(sib_matrix,1);
    
    frame = sib_matrix(i,4);
    first_cell = sib_matrix(i,2);
    second_cell = sib_matrix(i,3);
    
    bool1 = xyzs_id(:,xyzs_id_columns) == first_cell; 
    frame_vector1 = xyzs_id(bool1,xyzs_id_columns-1);
    sib_vector1 = xyzs_id(bool1,xyzs_id_columns-4);
    
    bool2 = xyzs_id(:,xyzs_id_columns) == second_cell;    
    frame_vector2 = xyzs_id(bool2,xyzs_id_columns-1);
    sib_vector2 = xyzs_id(bool2,xyzs_id_columns-4);

    % if the first frame the cell appears is the first frame in which it
    % has a sibling, the cell is a sibling (1); else assume the cell is a
    % collision (2) [next loop will correct for continuing collisions and
    % exits]
    if frame_vector1(1) == frame || frame_vector2(1) == frame
        sib_matrix(i,6) = 1;
    else
        sib_matrix(i,6) = 2;
    end
    
    new_bool1 = frame_vector1 < frame & frame_vector1 >= (frame-max_collision_time);
    new_bool2 = frame_vector2 < frame & frame_vector2 >= (frame-max_collision_time);
    
    new_sib_vector1 = sib_vector1(new_bool1);
    new_sib_vector2 = sib_vector2(new_bool2);
  
    if sum(new_sib_vector1) ~=0 && sum(new_sib_vector2) ~= 0
       new_index1 = find((sib_matrix(:,2) == first_cell | sib_matrix(:,3) == first_cell) & sib_matrix(:,4) < frame,1,'last'); 
       new_index2 = find((sib_matrix(:,2) == second_cell | sib_matrix(:,3) == second_cell) & sib_matrix(:,4) < frame,1,'last');
       
       if new_index1 == new_index2
           sib_matrix(i,6) = 3;
           sib_matrix(i,7) = new_index1;
       end    
    else
        % If there are siblings for either cell during the time a collision can
        % occur before the current sibling event, the event is listed as a
        % continuing collision
        if sum(new_sib_vector1) ~=0
           new_index = find((sib_matrix(:,2) == first_cell | sib_matrix(:,3) == first_cell) & sib_matrix(:,4) < frame,1,'last');
           sib_matrix(i,7) = new_index;
           second_collision_frame = sib_matrix(new_index,4);
           if sib_matrix(new_index,2) == second_cell | sib_matrix(new_index,3) == second_cell
               sib_matrix(i,6) = 3;
           end
        end

        if sum(new_sib_vector2) ~= 0
            new_index = find((sib_matrix(:,2) == second_cell | sib_matrix(:,3) == second_cell) & sib_matrix(:,4) < frame,1,'last');
            sib_matrix(i,7) = new_index;
            second_collision_frame = sib_matrix(new_index,4);
           if sib_matrix(new_index,2) == first_cell || sib_matrix(new_index,3) == first_cell
               sib_matrix(i,6) = 3;
           end
        end        
    end
    
%%%%% Code that was added to check if divisions are truly divisions, or are
%%%%% part of a previous merging event that is mislabeled    
        if frame_vector1(1) == frame || frame_vector2(1) == frame
            if frame_vector1(1) == frame % If the first cell is the newly dividing cell, then check the second cell's history
                div_bool2 = frame_vector2 < frame & frame_vector2 >= (frame-3*max_collision_time);
                new_sib_vector2 = sib_vector2(div_bool2);

                if sum(new_sib_vector2) ~=0 % If the second cell had a collision event, let's check to see if the division is really part of that previous event
                    new_index = find((sib_matrix(:,2) == second_cell | sib_matrix(:,3) == second_cell) & sib_matrix(:,4) < frame,1,'last'); % Find the previous event
                    if sib_matrix(new_index,2) == second_cell; % Determine which cell it had an event with and list that as "old_cell"
                        old_cell = sib_matrix(new_index,3);
                        old_cell_frame = sib_matrix(new_index,4);
                    elseif sib_matrix(new_index,3) == second_cell;
                        old_cell = sib_matrix(new_index,2);
                        old_cell_frame = sib_matrix(new_index,4);
                    end
                    bool_old_cell = xyzs_id(:,xyzs_id_columns) == old_cell & xyzs_id(:,xyzs_id_columns-1) > old_cell_frame; %& xyzs_id(:,xyzs_id_columns-1) <= frame + max_collision_time;
                    if sum(bool_old_cell) == 0 % If the old_cell has disappeared, then this is part of that previous event, so label the event type a "3"
                        sib_matrix(i,6) = 3;
                        sib_matrix(i,7) = new_index;
                        % Relabel cell in sib_matrix and xyzs_id
                        sib_matrix(sib_matrix(:,2)==first_cell,2)=old_cell;
                        sib_matrix(sib_matrix(:,3)==first_cell,3)=old_cell;
                        xyzs_id(xyzs_id(:,xyzs_id_columns)==first_cell,xyzs_id_columns) = old_cell;
                        wrong_division_count = wrong_division_count + 1;
                        wrong_division(wrong_division_count,1) = second_cell;
                        wrong_division(wrong_division_count,2) = old_cell;
                        wrong_division(wrong_division_count,3) = first_cell;
                        wrong_division(wrong_division_count,4) = old_cell_frame;
                    end
                end

            elseif frame_vector2(1) == frame
                div_bool1 = frame_vector1 < frame & frame_vector1 >= (frame-3*max_collision_time);
                new_sib_vector1 = sib_vector1(div_bool1);

                if sum(new_sib_vector1) ~=0 % If the first cell had a collision event, let's check to see if the division is really part of that previous event
                    new_index = find((sib_matrix(:,2) == first_cell | sib_matrix(:,3) == first_cell) & sib_matrix(:,4) < frame,1,'last');
                    if sib_matrix(new_index,2) == first_cell;% Determine which cell it had an event with and list that as "old_cell"
                        old_cell = sib_matrix(new_index,3);
                        old_cell_frame = sib_matrix(new_index,4);
                    elseif sib_matrix(new_index,3) == first_cell;
                        old_cell = sib_matrix(new_index,2);
                        old_cell_frame = sib_matrix(new_index,4);
                    end
                    bool_old_cell = xyzs_id(:,xyzs_id_columns) == old_cell & xyzs_id(:,xyzs_id_columns-1) > old_cell_frame; %& xyzs_id(:,xyzs_id_columns-1) <= frame + max_collision_time;
                    if sum(bool_old_cell) == 0 % If the old_cell has disappeared, then this is part of that previous event, so label the event type a "3"
                        sib_matrix(i,6) = 3;
                        sib_matrix(i,7) = new_index;
                        % Now relabel the cell ids in xyzs_id and
                        % sib_matrix
                        sib_matrix(sib_matrix(:,2)==second_cell,2)=old_cell;
                        sib_matrix(sib_matrix(:,3)==second_cell,3)=old_cell;
                        xyzs_id(xyzs_id(:,xyzs_id_columns)==second_cell,xyzs_id_columns) = old_cell;
                        wrong_division_count = wrong_division_count + 1; 
                        wrong_division(wrong_division_count,1) = first_cell;
                        wrong_division(wrong_division_count,2) = old_cell;
                        wrong_division(wrong_division_count,3) = second_cell;
                        wrong_division(wrong_division_count,4) = old_cell_frame;
                    end
                end            
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
end

% Determine the number of events
event_counter = 1;
for x = 1:size(sib_matrix,1);
    if sib_matrix(x,6) < 3
        sib_matrix(x,8) = event_counter;
        event_counter = event_counter + 1;
    else
        sib_matrix(x,8) = sib_matrix(sib_matrix(x,7),8);
    end 
end    
         

% Organize the sib_matrix data into a cell array with each cell containing
% the sibling information for a single overall event
sib_matrix = sortrows(sib_matrix,8);    
event_array = cell(max(sib_matrix(:,8)),1); 

for s = 1:max(sib_matrix(:,8));
    bool3 = sib_matrix(:,8) == s;
    event_array{s} = sib_matrix(bool3,1:9);
    
    % loop that goes from 2 to the number of frames for the given event
    % (last frame minus first frame)
    if size(event_array{s},1) > 1;  % if there is more than 1 frame      
        
        frame_length = event_array{s}(size(event_array{s},1),4) - event_array{s}(1,4);
        for t = 2:frame_length;
         
            if event_array{s}(t,4) ~= (event_array{s}(t-1,4) + 1); % if a frame is missing
                event_array{s} = insertrows(event_array{s},zeros(1,9),t - 1); % insert row of 0's after missing frame
                
                current_frame = event_array{s}(t-1,4) + 1; % update the current frame (missing frame)
                cell_index = find(event_array{s}(1:t,6) <= 4 & event_array{s}(1:t,6) >= 1, 1, 'last'); % index for first cell based on previous sibling event
                first_cell = event_array{s}(cell_index,2); % identify the first cell
                second_cell = event_array{s}(cell_index,3); % identifiy the second cell
                
                first_index = find(xyzs_id(:,xyzs_id_columns) == first_cell & xyzs_id(:,xyzs_id_columns-1) == current_frame);
                second_index = find(xyzs_id(:,xyzs_id_columns) == second_cell & xyzs_id(:,xyzs_id_columns-1) == current_frame);
                
                event_array{s}(t,4) = current_frame; % update the current frame into event_array
                event_array{s}(t,8) = s; % update the event number into event_array
                
                if isempty(first_index) == 0 && isempty(second_index) == 0; % both cells appear in the missing frame
                    event_array{s}(t,2) = first_cell;
                    event_array{s}(t,3) = second_cell;
                    event_array{s}(t,6) = 5;
                    event_array{s}(t,1) = find(xyzs_id(:,xyzs_id_columns)==first_cell & xyzs_id(:,xyzs_id_columns-1) == current_frame);
                    event_array{s}(t,9) = find(xyzs_id(:,xyzs_id_columns)==second_cell & xyzs_id(:,xyzs_id_columns-1) == current_frame);
                elseif isempty(first_index) == 0 && isempty(second_index) == 1; % first cell shows up in frame but not second cell
                    event_array{s}(t,2:3) = first_cell;
                    event_array{s}(t,6) = 6;
                elseif isempty(first_index) == 1 && isempty(second_index) == 0; % second cell shows up in frame but not first cell
                    event_array{s}(t,2:3) = second_cell;
                    event_array{s}(t,6) = 6;
                elseif isempty(first_index) == 1 && isempty(second_index) == 1; % neither cell shows up in current frame
                    event_array{s}(t,6) = 7;
                end
            end
        end
    end                  
end