function [xyzs_id event_array missing_cell_mat] = mult_division_adjuster(xyzs_id,xyzs_id_columns,event_array,mult_event_array, max_collision_time)

% FUNCTION 1: Check 2-body merging events and see if they stemmed from
% complex merging event. If they did, then make sure no cells in that
% complex event disappeared
min_frame = max_collision_time;

missing_index = 1; %Initialize the missing cell index
missing_cell_mat = zeros(size(event_array,1),3);
count = 0;
count1 = 0;
count2 = 0;

for i=1:size(event_array,1)
    % Check if the event is a division and occurs after the min number of
    % frames to consider a division
    if event_array{i}(1,6) == 1 && event_array{i}(1,4) >= min_frame
        count2 = count2 + 1;
        % Find out which cell is the newly dividing cell
        first_cell = event_array{i}(1,2);
        first_cell_xyzs = xyzs_id(xyzs_id(:,xyzs_id_columns)==first_cell,:);
        second_cell = event_array{i}(1,3);
        second_cell_xyzs = xyzs_id(xyzs_id(:,xyzs_id_columns)==second_cell,:);
        division_frame = event_array{i}(1,4);
        if first_cell_xyzs(1,xyzs_id_columns - 1) == division_frame; %If first cell is new cell
            new_cell = first_cell;
            old_cell = second_cell;
            old_cell_xyzs = second_cell_xyzs;
        elseif second_cell_xyzs(1,xyzs_id_columns - 1) == division_frame; %If second cell is new cell
            new_cell = second_cell;
            old_cell = first_cell;
            old_cell_xyzs = first_cell_xyzs;
        end
        
        % Check to see if old cell had a complex merging event within
        % max_collision_time of the division
        complex_frame_bool = old_cell_xyzs(:,xyzs_id_columns-1) >= (division_frame-max_collision_time) & old_cell_xyzs(:,xyzs_id_columns-1) < division_frame;
        complex_check = old_cell_xyzs(complex_frame_bool,8); % Note column 8 is the complex tag
        % If it did have a complex event, find that event and make sure
        % no cells were lost
        if sum(complex_check) > 0
            count1 = count1 + 1;
            % Find potential events it is involved in based on frame
            event_end_frame_vec = cell2mat(mult_event_array(:,4));
            bool_event = event_end_frame_vec < division_frame & event_end_frame_vec >= (division_frame - max_collision_time);
            event_index = [1:length(event_end_frame_vec)];
            event_index = event_index(bool_event);
            % Take potential events and find which event it is involved in
            for j = 1:length(event_index)
                cell_list = mult_event_array{event_index(j),2};
                % If it is part of this event, check other cells to see if
                % they disapper
                if ismember(old_cell,cell_list)
                    for k=1:length(cell_list)
                        cell_xyzs = xyzs_id((xyzs_id(:,xyzs_id_columns) == cell_list(k)&xyzs_id(:,xyzs_id_columns-1) > division_frame),:);
                        % If the cell has disappeared, reclassify this
                        % event as "11" instead of "1", highlighting it is
                        % likely not a division
                        if isempty(cell_xyzs)
                            event_array{i,6} = 11;
                            count = count + 1;
                            missing_cell = cell_list(k);
                            missing_cell_mat(missing_index,1) = missing_cell; % Lost cell ID
                            missing_cell_mat(missing_index,2) = new_cell; % Incorrectly labeled division cell
                            missing_cell_mat(missing_index,3) = division_frame; % Frame this relabeling occurred
                            missing_cell_mat(missing_index,4) = i;
                            missing_index = missing_index + 1;
                        end
                    end
                end                
            end
        end
    end
end

missing_cell_mat = missing_cell_mat(missing_cell_mat(:,1)~=0,:);

% FUNCTION 2: Check multi-body merging events for new cells that appear.
% Then check to see if cells in that event have disappeared.  If not, then
% classify this as a division.


