function [xyzs_id, mastertable] = collision_tags_rtc2_3(xyzs_id, event_array, xyzs_id_columns, frame_avg)
%   function [xyzs_id, sib_lookup, cost_info, ecount, straight_count, crossover_count, sib_cost_info, cost_array, mastertable] = collision_tags_rtc2_3(xyzs_id, event_array, xyzs_id_columns, frame_avg)
%
%   This function 'fixes' mislabeled ID tags resulting from cell merging
%   events. It achieves this by running a positional cost function on the
%   cells involved in a given collision and calculates the cost for
%   switching cell IDs after the collision. 
%
%  2/1/2013
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
%  frame_avg: matrix of particle characteristic averages for each frame
%
%  OUTPUTS:
%  xyzs_id: updated matrix of particle information after 
%  mastertable: table used for relabeling cells after doing merging event
%           analysis
%
%  LOCAL PARAMETERS:

%Check whether collision corrector has been run - if not, run it
if nargin < 2
    error('Please run collision corrector before continuing.');
end

%Note: The variable 'event_array' contains all of the particle information 
%for sorting out collisions
nume = size(event_array,1); %Number of events
xyzs_id(:,xyzs_id_columns+3) = zeros(size(xyzs_id,1),1);   %Insert a new column for retagging
sib_lookup = zeros(2*nume,3);   %Create a look-up table based on cost function
scount = 1;
cost_info = zeros(2*nume,8);
cost_count = 1; 
ecount = 1;
sib_cost_info = zeros(nume,6);
% sib_cost_info contains: 1) event_number, 2) cell 1 involved, 3) cell 2 
% involved, 4) straight/crossover identified by cost function, 5) end cell 
% id 1, 6) end cell id 2 
straight_count = 0; %Count number of straight line cases classified by cost function. 
crossover_count = 0; %Count number of crossover cases identified by cost function. 

cost_array = cell(nume,1);
mastertable = zeros(nume,9);
event_in_pos_vec = zeros(nume,1);
final_frame = max(xyzs_id(:,xyzs_id_columns-1));
no_cost_counter = 0;


for i = 1:nume
    event_info = event_array{i};    %Event information
    esize = size(event_info,1);     %Size of event info
    p1 = event_info(1,2);           %Cell 1 ID
    p2 = event_info(1,3);           %Cell 2 ID
    frame_start = event_info(1,4);  %Frame that event starts in
    ctype = event_info(1,6);
    sib_cost_info(i,1) = i;
    sib_cost_info(i,2) = p1;
    sib_cost_info(i,3) = p2;
    
    %If the event occurs in the final frame skip all the cost analysis
    if frame_start == final_frame
        no_cost_counter = no_cost_counter + 1;
        no_cost_vector(no_cost_counter,1) = i;
        continue
    end
    
    %Check if entry is only one row in size:
    if esize == 1
        ecount = ecount+1;  %Counter for number of single entry events
        
        %Check for division case
        if ctype == 1
            prow_index = (xyzs_id(:,xyzs_id_columns) == event_info(1,2)) & (xyzs_id(:,xyzs_id_columns-1) == event_info(1,4));
            srow_index = (xyzs_id(:,xyzs_id_columns) == event_info(1,3)) & (xyzs_id(:,xyzs_id_columns-1) == event_info(1,4));
            xyzs_id(prow_index,xyzs_id_columns+3) = p1;
            xyzs_id(srow_index,xyzs_id_columns+3) = p2;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This is put in to check events/cost_info for 
            % debugging/optimization 
            %
              cost_info(cost_count,:) = 1;
              cost_count = cost_count+2;
              sib_lookup(scount,:) = 1;
              sib_lookup(scount+1,:) = 1;
              scount = scount+1; 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            continue
        end
        
        %Cell information used for the "entering" cells will be obtained
        %from the frame of first sibling identification
        xyzs_id_c1 = xyzs_id(event_info(1,1),:);
        xyzs_id_c2 = xyzs_id(event_info(1,9),:);
        
        %Cell information for the "exiting" cell will be the first frame 
        %after the collision event where both cells are present; this finds
        %the first frame after collision where they both are present
        bool_cell1 = xyzs_id(:,xyzs_id_columns) == p1 & xyzs_id(:,xyzs_id_columns - 1) > frame_start;
        bool_cell2 = xyzs_id(:,xyzs_id_columns) == p2 & xyzs_id(:,xyzs_id_columns - 1) > frame_start;
        cell1_info = xyzs_id(bool_cell1,:);
        cell2_info = xyzs_id(bool_cell2,:);
        
        v=1; %initialize counter for loop
        exit_frame = []; %create empty exit_frame matrix to enter while loop
        while isempty(exit_frame) && v <= size(cell1_info,1) %loop finds frame where both cells are first present
            exit_frame_index = find(cell2_info(:,xyzs_id_columns-1) == cell1_info(v,xyzs_id_columns-1));
            exit_frame = cell2_info(exit_frame_index,xyzs_id_columns-1);
            v = v+1;
        end
        
        if isempty(exit_frame)
            no_cost_counter = no_cost_counter + 1;
            no_cost_vector(no_cost_counter,1) = i;
            continue
        end
        
        bool_cell1_exit = cell1_info(:,xyzs_id_columns-1) == exit_frame;
        xyzs_id_c3 = cell1_info(bool_cell1_exit,:);
        xyzs_id_c4 = cell2_info(exit_frame_index,:);
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %Check that there is particle information after event frame:
        if isempty(xyzs_id_c1) || isempty(xyzs_id_c2) || isempty(xyzs_id_c3) || isempty(xyzs_id_c4)
            %Fill in corresponding rows in cost_info and sib_lookup with
            %1's for now --> this will have to change. Placeholder to check
            %some results from cost function. 
            cost_info(cost_count,:) = 2;
            cost_count = cost_count+2;
            sib_lookup(scount,:) = 2;
            sib_lookup(scount+1,:) = 2;
            scount = scount+2; 
            continue
        end
        
            % Calculate the distance between centroids
        cost13 = ((xyzs_id_c1(1) - xyzs_id_c3(1))^2 + (xyzs_id_c1(2) - xyzs_id_c3(2))^2)^(1/2); 
        cost24 = ((xyzs_id_c2(1) - xyzs_id_c4(1))^2 + (xyzs_id_c2(2) - xyzs_id_c4(2))^2)^(1/2);
        cost1 = cost13 + cost24;

        cost14 = ((xyzs_id_c1(1) - xyzs_id_c4(1))^2 + (xyzs_id_c1(2) - xyzs_id_c4(2))^2)^(1/2); 
        cost23 = ((xyzs_id_c2(1) - xyzs_id_c3(1))^2 + (xyzs_id_c2(2) - xyzs_id_c3(2))^2)^(1/2);
        cost2 = cost14+cost23;
        
        cell1_start_index = event_info(1);
        cell2_start_index = event_info(9);
        cell1_end_index = find(xyzs_id(:,xyzs_id_columns-1) == exit_frame & xyzs_id(:,xyzs_id_columns) == p1);
        cell2_end_index = find(xyzs_id(:,xyzs_id_columns-1) == exit_frame & xyzs_id(:,xyzs_id_columns) == p2);
        
        %Make master lookup table
        mastertable(i,1)=xyzs_id_c1(xyzs_id_columns-1);
        mastertable(i,2)=xyzs_id_c3(xyzs_id_columns-1);
        mastertable(i,3)=xyzs_id_c1(xyzs_id_columns);
        mastertable(i,4)=xyzs_id_c2(xyzs_id_columns);
        mastertable(i,5)=cell1_start_index;
        mastertable(i,6)=cell2_start_index;
        mastertable(i,7)=cell1_end_index;
        mastertable(i,8)=cell2_end_index;
        if cost1<cost2
            mastertable(i,9)=0;  % IDs are correct, don't switch
        else
            mastertable(i,9)=1; % Need to switch IDs
        end
        continue
    end
    
    %Loop through each entry and re-tag individually (new tag in xyzs_id)
    for j=1:esize
        if ctype == 1
            prow_index = (xyzs_id(:,xyzs_id_columns) == event_info(1,2)) & (xyzs_id(:,xyzs_id_columns-1) == event_info(1,4));
            srow_index = (xyzs_id(:,xyzs_id_columns) == event_info(1,3)) & (xyzs_id(:,xyzs_id_columns-1) == event_info(1,4));
            xyzs_id(prow_index,xyzs_id_columns+3) = p1;
            xyzs_id(srow_index,xyzs_id_columns+3) = p2;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %This will also have to be replaced. Temporarily put in to
            %check events/cost_info 
            
            cost_info(cost_count,:) = 1;
            cost_count = cost_count+2;
            sib_lookup(scount,:) = 1;
            sib_lookup(scount+1,:) = 1;
            scount = scount+1; 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            continue
        end
        
            if event_info(j,6) ~= 7
            prow_index = (xyzs_id(:,xyzs_id_columns) == event_info(j,2)) & (xyzs_id(:,xyzs_id_columns-1) == event_info(j,4));
            srow_index = (xyzs_id(:,xyzs_id_columns) == event_info(j,3)) & (xyzs_id(:,xyzs_id_columns-1) == event_info(j,4));
            xyzs_id(prow_index,xyzs_id_columns+3) = p1;
            xyzs_id(srow_index,xyzs_id_columns+3) = p2;

%update for cost function
            %Check whether the ending particle is the same as the beginning
            %particle - retag if not
            if j == esize                
                %Use position cost function for all cases 
                event_in_pos_vec(i) = 1;
                [ lookup_table, pos_matrix ] = pos_2( event_info, xyzs_id, xyzs_id_columns, frame_avg );             

                    mastertable(i,1) = lookup_table(1,1);
                    mastertable(i,2) = lookup_table(3,1);
                    mastertable(i,3) = lookup_table(1,2);
                    mastertable(i,4) = lookup_table(2,2);
                    mastertable(i,5) = lookup_table(1,4);
                    mastertable(i,6) = lookup_table(2,4);
                    mastertable(i,7) = lookup_table(3,4);
                    mastertable(i,8) = lookup_table(4,4);
                    if lookup_table(3,2) == lookup_table(3,3)
                        mastertable(i,9) = 0;
                    else
                        mastertable(i,9) = 1;
                    end
            end
        end
    end
end


mastertable(all(mastertable==0,2),:)=[];
mastertable = sortrows(mastertable, -2);

end
