function [xyzs_id event_array div_angle n_divisions] = division_corrector(xyzs_id, event_array, xyzs_id_columns, particle_radius, image, div_toggle, max_collision_time, border_width, frust_toggle, stackname)
% Function to analyze interaction events classified as divisions and
% reclassify incorrect divisions. The code then duplicates the parent
% information of a dividing cell so both daughter cells have the parent
% information prior to the division.
%
%  6/05/2013
%  R. Baker, M. Brasch
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
%  xyzs_id: matrix of particle information post-tracking 
%  event_array: cell_array output from collision_corrector.m that contains
%               information for every event
%  xyzs_id_columns: column of the xyzs_id matrix containing cell ID info
%  particle_radius: radius of the average cell
%  image: image matrix (used for size of the image)
%  min_frame: first frame to start counting divisions
%  border_width: number of pixels away from image borders to not count 
%                divisions
%  frust_toggle: 1-on; 2-off  (frustrated divisions that don't completely
%                divide and may cause many more false positive divisions 
%                are corrected)
%
%  OUTPUTS:
%  xyzs_id: matrix of particle information after duplicating parent
%           information of divisions
%  event_array: updated event_array after correcting for frustrated
%               divisions and reclassifying incorrect division events
%
%  LOCAL PARAMETERS:

% Number of frames after a division to see if cells are far enough away to
% have fully divided (frustrated division analysis)
frust_frame_min = 40;
tolerance = 0.7; % Difference area average intensity can be for cells to still be considered dividing cells
min_frame = max_collision_time;

symmetric_div_toggle = 1;
% border_width = 30; %%%%%%%%%%%%%%%%%%%%%
% particle_radius = 10;%%%%%%%%%%%%%%%%%%%%%%
% image = A;%%%%%%%%%%%%%%%%%%%%%

if nargin < 9
    frust_toggle = 0;
end

if nargin < 8
    border_width = 30; % distance away from border an interaction event must be to be considered a division
end

if nargin < 7
    min_frame = 10; % first frame an interaction event can be considered a division
end


% Generate a division matrix that contains (columns): 1-cell 1, 2-cell 2, 
% 3-frame number, 4-event number, 5-Is it a frustrated division
div_matrix = zeros(size(event_array,1),5);
index = 1;
count1 = 0;
count2 = 0;
event2 = zeros(200,1);
count3 = 0;
event3 = zeros(200,1);
% When reclassifying divisions: 
   %Check 1 (change to 11) = occurs before min frame 
   %Check 2 (change to 12) = too close to border
   %Check 3 (change to 13) = does not meet area criterion (only if
   %                         symmetric division)
for i=1:size(event_array,1);
    event_info = event_array{i};
    if event_array{i}(1,6) == 1 % If it is originally classified as an event
        % Extract useful information to be used for the next division
        % checks
        cell1_bool = xyzs_id(:,xyzs_id_columns)==event_info(1,2) & xyzs_id(:,xyzs_id_columns-1)==event_info(1,4);
        cell2_bool = xyzs_id(:,xyzs_id_columns)==event_info(1,3) & xyzs_id(:,xyzs_id_columns-1)==event_info(1,4);
        cell1_mat = xyzs_id(cell1_bool,:);
        cell2_mat = xyzs_id(cell2_bool,:);
        
        
        %%%%%%%%%%%%%%%%%%  CHECK 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check the criterion of occuring after min_frame
        if event_info(1,4) < min_frame; 
            event_array{i}(1,6) = 11; % reclassify division event as before min frame
            count1 = count1+1; % number of ignored division events
            continue % move to next division event
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%  CHECK 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check the criterion of occuring within the border
        pos_mat = [cell1_mat(1:2); cell2_mat(1:2)];
        pos_bool = pos_mat(:,1) < border_width | pos_mat(:,2) < border_width | pos_mat(:,2) > (size(image,1)-border_width) | pos_mat(:,1) > (size(image,2)-border_width);
        if sum(pos_bool) > 0
            event_array{i}(1,6) = 12; % reclassify division event as outside border
            count2 = count2+1; % number of ignored division events
            event2(count2) = i;
            continue % move to next division event
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%  CHECK 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check the criterion of similar area (only if symmetric division)
        if symmetric_div_toggle == 1;
            area1 = cell1_mat(6);
            area2 = cell2_mat(6);
            if abs(area1-area2)/mean([area1 area2]) > tolerance
                [tolerance_check sibling_cell] = div_param_check(xyzs_id, xyzs_id_columns, event_info, cell2_mat, tolerance, 4*particle_radius);
                if tolerance_check == 0;
                    event_array{i}(1,6) = 13; % reclassify division event as failing area criterion
                    count3 = count3+1; % number of ignored division events
                    event3(count3) = i;
                    continue % move to next division event
                else
                    event_array{i}(1,2) = sibling_cell;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        div_matrix(index,1) = event_array{i}(1,2);
        div_matrix(index,2) = event_array{i}(1,3);
        div_matrix(index,3) = event_array{i}(1,4);
        div_matrix(index,4) = i;
        div_matrix(index,5) = 0;
        
        index = index + 1;
    end
end
div_matrix(div_matrix(:,1)==0,:) = [];

if ~isempty(div_matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   frustrated division analysis
    if frust_toggle == 1;
        % Check to see if these are frustrated divisions based on particle distance
        % away from each other, or if either cell was previously involved in a
        % prior frustrated division
        div_matrix = sortrows(div_matrix,3);
        dist_vec = zeros(size(div_matrix,1),3);
        count = 0;
        index = 1;
        for j=1:size(div_matrix,1);
            % Find cells that have been classified as frustrated divisions
            % (which is determined later in this loop)
            repeat_frust_bool = div_matrix(:,1) == div_matrix(j,1) | div_matrix(:,1) == div_matrix(j,2) | div_matrix(:,2) == div_matrix(j,1) | div_matrix(:,2) == div_matrix(j,2);
            if sum(div_matrix(repeat_frust_bool,5)) > 0 % If either cell had a previous frustrated division, count this as a continuation of that division
                div_matrix(j,5) = 2; % continuation of a frstrated division
                count = count + 1;
            else % Check to see if this is a frustrated division
                cell1_bool = xyzs_id(:,xyzs_id_columns) == div_matrix(j,1) & xyzs_id(:,xyzs_id_columns-1) > div_matrix(j,3);
                cell2_bool = xyzs_id(:,xyzs_id_columns) == div_matrix(j,2) & xyzs_id(:,xyzs_id_columns-1) > div_matrix(j,3);
                cell1_mat = xyzs_id(cell1_bool,:);
                cell2_mat = xyzs_id(cell2_bool,:);
                % find the frames the two cells in the division share
                [shared_frames, index1, index2] = intersect(cell1_mat(:,xyzs_id_columns-1), cell2_mat(:,xyzs_id_columns-1));
                if isempty(shared_frames)
                    check_event(index) = div_matrix(j,4);
                    index = index+1;
                    continue
                end

                if size(shared_frames,1) > frust_frame_min
                    frame_calc = frust_frame_min;
                else
                    frame_calc = size(shared_frames,1);
                end

                % Calculate distance between the two cells for every frame they
                % share
                dist = zeros(frame_calc,1);            
                for r=1:frame_calc
                    dist(r) = sqrt((cell1_mat(index1(r),1) - cell2_mat(index2(r),1))^2 + (cell1_mat(index1(r),2) - cell2_mat(index2(r),2))^2);
                end

                dist_vec(j,1) = max(dist);
                dist_vec(j,2) = div_matrix(j,4);
                dist_vec(j,3) = r;

                % If the cells never move away, it is a frustrated division
                if max(dist) <= 2*particle_radius
                    div_matrix(j,5) = 1;
                end
            end
        end

        dist_vec(dist_vec(:,2)==0,:)=[];

        frust_bool = div_matrix(:,5) == 1;
        frust_matrix = div_matrix(frust_bool,:);

        double_frust_bool = div_matrix(:,5) == 2;
        double_frust_matrix = div_matrix(double_frust_bool,:);
        div_switch_matrix = zeros(size(double_frust_matrix,1),3);

        % Generate a relabeling matrix to get rid of new particles attributed to
        % multiple frustrated divisons
        for x = 1:size(double_frust_matrix,1);
            cell1 = double_frust_matrix(x,1);
            cell2 = double_frust_matrix(x,2);
            frame = double_frust_matrix(x,3);

            prior_bool = div_matrix(:,1) == cell1 | div_matrix(:,2) == cell1 | div_matrix(:,1) == cell2 | div_matrix(:,2) == cell2 & div_matrix(:,3) < frame;
            prior_matrix = div_matrix(prior_bool,:);

            if cell1 == prior_matrix(1,1) | cell1 == prior_matrix(1,2);
                if cell1 == prior_matrix(1,1);
                    new_cellID = prior_matrix(1,2);
                    old_cellID = cell2;
                else
                    new_cellID = prior_matrix(1,1);
                    old_cellID = cell2;
                end

                % This makes sure that the ID we are switching to doesn't already
                % exist, which it may if the frustrated division divides and then
                % has another frustrated division, or if the division is due to
                % some complex merging event mistake.
                bool_check = xyzs_id(:,xyzs_id_columns) == new_cellID & xyzs_id(:,xyzs_id_columns-1) > frame;
                if sum(bool_check) ==0;
                    div_switch_matrix(x,1) = old_cellID; % old ID
                    div_switch_matrix(x,2) = new_cellID; % ID to switch to
                    div_switch_matrix(x,3) = frame; % frame to begin switch
                    div_matrix(div_matrix(:,1) == old_cellID, 1) = new_cellID; % update div matrix with new IDs
                    div_matrix(div_matrix(:,2) == old_cellID, 2) = new_cellID;
                    double_frust_matrix(double_frust_matrix(:,1) == old_cellID, 1) = new_cellID;
                    double_frust_matrix(double_frust_matrix(:,2) == old_cellID, 2) = new_cellID;
                end

            elseif cell2 == prior_matrix(1,1) | cell2 == prior_matrix(1,2);
                if cell2 == prior_matrix(1,1);
                    new_cellID = prior_matrix(1,2);
                    old_cellID = cell1;
                else
                    new_cellID = prior_matrix(1,1);
                    old_cellID = cell1;
                end

                bool_check = xyzs_id(:,xyzs_id_columns) == new_cellID & xyzs_id(:,xyzs_id_columns-1) > frame;
                if sum(bool_check) == 0;
                    div_switch_matrix(x,1) = old_cellID; % old ID
                    div_switch_matrix(x,2) = new_cellID; % ID to switch to
                    div_switch_matrix(x,3) = frame; % frame to begin switch
                    div_matrix(div_matrix(:,1) == old_cellID, 1) = new_cellID; % update div matrix with new IDs
                    div_matrix(div_matrix(:,2) == old_cellID, 2) = new_cellID;
                    double_frust_matrix(double_frust_matrix(:,1) == old_cellID, 1) = new_cellID;
                    double_frust_matrix(double_frust_matrix(:,2) == old_cellID, 2) = new_cellID;
                end
            end
        end

        div_switch_matrix(div_switch_matrix(:,1)==0,:) = [];
        
        % Let's correct the frustrated cell IDs in the xyzs_id matrix (only if
        % you have true frustrated divisions)
        if frust_toggle ==1;
            for j = 1:size(div_switch_matrix,1);
                xyzs_id(xyzs_id(:,xyzs_id_columns) == div_switch_matrix(j,1) & xyzs_id(:,xyzs_id_columns-1) >= div_switch_matrix(j,3),xyzs_id_columns) = div_switch_matrix(j,2);
            end
        end
        
        % Update event_array for ids that may have changed
        for s=1:size(event_array,1);
            event_info = event_array{s};
            unique1 = unique(event_info(:,2));
            unique2 = unique(event_info(:,3));
            unique_cells = unique([unique1; unique2]);

            for u = 1:length(unique_cells);
                bool = div_switch_matrix(:,1) == unique_cells(u);
                if sum(bool) > 0
                    old_cellID = unique_cells(u);
                    new_cellID = div_switch_matrix(bool,2);
                    bool_event_info1 = event_info(:,2) == old_cellID;
                    bool_event_info2 = event_info(:,3) == old_cellID;
                    if sum(bool_event_info1) > 0;
                        event_info(bool_event_info1,2) = new_cellID;
                    end
                    if sum(bool_event_info2) > 0;
                        event_info(bool_event_info2,3) = new_cellID;
                    end
                end    
            end
            event_array{s} = event_info;
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initiate division GUI
if div_toggle == 1    
   [div_output, pause_status, div_event] = division_gui(stackname, event_array, xyzs_id, xyzs_id_columns);   
   
   %Check to make sure the user did not pause the GUI - if they did, ignore
   %output
   if pause_status == 0
       %Incorporate GUI output into division matrix
       final_div_output = [div_output, div_event];
       bool_zeros = final_div_output(:,1) == 0;
       event_change = final_div_output(bool_zeros,2);

       if ~isempty(event_change)
            for i=1:size(event_change,1)
                cell_id = event_change(i);
                div_matrix(div_matrix(:,4) == cell_id,:) = 0;
            end
            div_matrix(div_matrix(:,1)==0,:) = [];
       end
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Now we want to run the duplicate code for division analysis:
    if frust_toggle == 1;
        div_dup_bool = div_matrix(:,5) ~= 2; % don't count multiple frustrated divisions
    else
        div_dup_bool = ones(div_matrix,1); % count all divisions
    end
    
    div_dup_matrix = div_matrix(div_dup_bool,:);
    count = 0;
    weird = 0;
    
    div_angle = zeros(size(div_dup_matrix,1),1);
    index1 = 1;
    % Duplicate the parent information and add it to xyzs_id
    for y=1:size(div_dup_matrix,1);
        cell1 = div_dup_matrix(y,1); % first cell in division
        cell2 = div_dup_matrix(y,2); % second cell in division
        frame = div_dup_matrix(y,3); % frame division begins
        
        cell1_bool = xyzs_id(:,xyzs_id_columns)==cell1 & xyzs_id(:,xyzs_id_columns-1)==frame;
        cell2_bool = xyzs_id(:,xyzs_id_columns)==cell2 & xyzs_id(:,xyzs_id_columns-1)==frame;
        cell1_mat = xyzs_id(cell1_bool,:);
        cell2_mat = xyzs_id(cell2_bool,:);
        if isempty(cell1_mat) | isempty(cell2_mat)
        else
            div_angle(index1,1) = div_dup_matrix(y,4);
            div_angle(index1,2) = atan2((cell1_mat(1,2)-cell2_mat(1,2)),(cell1_mat(1,1)-cell2_mat(1,1)));
            index1 = index1+1;
        end

        xyzs_bool1 = xyzs_id(:,xyzs_id_columns) == cell1 & xyzs_id(:,xyzs_id_columns-1)<frame;
        xyzs_bool2 = xyzs_id(:,xyzs_id_columns) == cell2 & xyzs_id(:,xyzs_id_columns-1)<frame;

        % Determine which cell daughter cell already has the parent
        % information prior to the division and duplicate and assign the
        % information to other daughter cell
        if sum(xyzs_bool1) > 0 && sum(xyzs_bool2) == 0
            xyzs_dup_mat = xyzs_id(xyzs_bool1,:);
            xyzs_dup_mat(:,xyzs_id_columns) = cell2;
            xyzs_id = vertcat(xyzs_id, xyzs_dup_mat);
        elseif sum(xyzs_bool2) > 0 && sum(xyzs_bool1) == 0
            xyzs_dup_mat = xyzs_id(xyzs_bool2,:);
            xyzs_dup_mat(:,xyzs_id_columns) = cell1;
            xyzs_id = vertcat(xyzs_id, xyzs_dup_mat);
        end
    end

    % Sortxyzs_id before output
    xyzs_id = sortrows(xyzs_id, [xyzs_id_columns, xyzs_id_columns-1]);
end

chk_div_angle = exist('div_angle','var');
if chk_div_angle == 1
    div_angle(div_angle(:,1)==0,:) = [];
    n_divisions = size(div_dup_matrix,1);
else 
    %no divisions
    div_angle = [];
    n_divisions = [];
end