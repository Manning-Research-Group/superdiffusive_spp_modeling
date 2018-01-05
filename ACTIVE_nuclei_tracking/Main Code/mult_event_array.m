function [multi_event_array] = mult_event_array(mult_sib_array,max_collision_time)

%mult_mat entries: 1) frame number, 2) event number, 3) flag (yes/no) -
%already processed event, 4) individual index number, 5) overall event
%number
mult_mat = [mult_sib_array{:,1};mult_sib_array{:,2};mult_sib_array{:,4}]';
mult_mat(:,4) = 1:size(mult_mat,1);
mult_mat(:,5) = zeros(size(mult_mat,1),1);

%Set up event array to hold multi-cell interaction information
multi_event_array = cell(size(mult_mat,1),4);

weight = 0.5;   %Number of cells that need to be present in a second mult-body interaction
event = 1;      %Initialize event number

%Main loop to sort out continuing events
for i = 1:size(mult_sib_array,1)
    % First check if this event has already been processed
    if mult_mat(i,3) == 0
        mult_mat(i,3) = 1;      %Update to indicate event is processed
        mult_mat(i,5) = event;  %Include event number
        cell_id = mult_sib_array{i,3};  %Find all cells involved in this event
        start_frame = mult_sib_array{i,1};  %Identify frame that this multi-body interaction occurs in
        check = 0;  % Check to identify continuing events
        while check == 0 
            %Check that the same cells do not have another multi-body merging event within max_collision time 
            bool_event = mult_mat(:,1)>start_frame & mult_mat(:,1) <= start_frame+max_collision_time;
            sib_bool_array = mult_mat(bool_event,:);
            check = 1;
            %Check all events for the same IDs - must meet threshold cell
            %number to be considered same event number
            for j = 1:size(sib_bool_array,1)
               temp_id = mult_sib_array{sib_bool_array(j,4),3};
               diff_id = setdiff(cell_id,temp_id);
               thresh = ceil(weight*length(cell_id));
               %Check if any new IDs have shown up which should now part be
               %part of the multi-body event
               if length(diff_id)<thresh
                   mult_mat(sib_bool_array(j,4),3) = 1;
                   mult_mat(sib_bool_array(j,4),5) = event;
                   new_id = setdiff(temp_id,cell_id);
                   if ~isempty(new_id)
                        cell_id = [cell_id;new_id];
                   end
                   start_frame = mult_mat(sib_bool_array(j,4),1);   %update frame number and go again
                   check = 0;   %reset check so while loop continues
                   break
               end
            end   
        end
        
        %Update overall array and advance event number
        multi_event_array{event,1} = event;
        multi_event_array{event,2} = cell_id;
        event = event+1;
    end
end

%Add start and end frame numbers into multi_event_array
mult_mat = sortrows(mult_mat,5);
for i = 1:max(mult_mat(:,5))
    bool = mult_mat(:,5) == i;
    multi_event_array{i,3} = min(mult_mat(bool,1));
    multi_event_array{i,4} = max(mult_mat(bool,1));
end

%Delete empty cell rows
multi_event_array( all(cellfun(@isempty,multi_event_array),2), : ) = [];
