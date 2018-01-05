function [ sib_matrix, xyzs_id ] = sib_matrix_creation( xyzs_id, xyzs_id_columns )
%  Creates a matrix of information for all sibling events to be used for
%  generating events, as well as organizing the sibling events for color
%  coding during the video plotting of ellipses
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
%  xyzs_id: matrix of particle information post-tracking (output from
%           trackmem_new.m) or post-correcting (output from
%           collision_tags.m)
%  xyzs_id_columns: column of xyzs_id matrix that contains the cell ID info
%
%  OUTPUTS:
%  sib_matrix: matrix of sibling information including color coding
%  xyzs_id: updated xyzs_id matrix incorporating color information from
%           sib_matrix for plotting
%
%  LOCAL PARAMETERS:

% Loop to extract sibling information
r=1; % sibling matrix counter    
sib_size = length(find(xyzs_id(:,xyzs_id_columns-3)~=0));
sib_matrix = zeros(sib_size, 5);

for i = 1:length(xyzs_id)

    if xyzs_id(i,xyzs_id_columns-3) ~= 0
        sib_index = find(xyzs_id(:,xyzs_id_columns-2) == xyzs_id(i,xyzs_id_columns-3)); %find where unique 
        %sibling ID is equal to individual particle ID number
        cell_id = xyzs_id(sib_index, xyzs_id_columns); %sibling particle number 
        frame_id = xyzs_id(sib_index, xyzs_id_columns-1); %sibling frame ID
        current_cell_id = xyzs_id(i,xyzs_id_columns); %current cell ID
        % organize information into sibling matrix:
        sib_matrix(r,1) = i;
        sib_matrix(r,2) = current_cell_id; %current cell
        sib_matrix(r,3)= cell_id; %sibling cell
        sib_matrix(r,4)= frame_id; %frame division occurs
        r=r+1;
    end

end
disp('Finished extracting sibling information');
save('variables.mat', 'xyzs_id', 'sib_matrix'); 

if isempty(sib_matrix) == 1 % Check if siblings exist
    disp('No siblings found')
    xyzs_id(:,xyzs_id_columns+1) = 0; % add column of zeros for color index
    xyzs_id(:,xyzs_id_columns+2) = 0; % add column of zeros for color index
else    %Sort sibling information
    disp('Sorting sibling information')

    %code for plotting siblings
    color = 1;
    sib_matrix(1,5) = color;
    % 1st column of house_keep = cell ID number (only includes siblings); 2nd
    % column = unique color for each sibling pair
    house_keep(1,1) = sib_matrix(1,2);
    house_keep(2,1) = sib_matrix(1,3);
    house_keep(1,2) = color;
    house_keep(2,2) = color;

    %assign a unique color (as a unique number) for each sibling pair
    for j = 2:size(sib_matrix,1)    %changed from length(sib_matrix) to size(...) 4/17/12

        % Determine if the current cell id has already been given a color based
        % on the previous cell ids; if yes the index is assigned to that of the
        % previous cell id; if no the index is empty
        index = find(house_keep(:,1) == sib_matrix(j,2), 1,'first');

        % For cell ids not given a previous color, advance the color and assign
        % the cell id and its sibling pair the new color.  Also add this into
        % the house_keep matrix


        if isempty(index)
            color=color+1; 
            sib_matrix(j,5) = color;

            house_keep(2*j-1,2) = color;
            house_keep(2*j,2) = color;

        % For cell ids that have already been given a color, set the sibling
        % color to the already assigned color
        else
            sib_matrix(j,5) = house_keep(index,2); %sibling color
            house_keep(2*j-1,2) = house_keep(index,2);
            house_keep(2*j,2) = house_keep(index,2);
        end    
        house_keep(2*j-1,1) = sib_matrix(j,2);
        house_keep(2*j,1) = sib_matrix(j,3);

    end

    save('variables.mat', 'xyzs_id', 'sib_matrix'); 


    % initialize the counter k and the color_index
    k=1;
    color_index = 1;

    %create a color matrix the length of the number of unique colors (unique
    %sibling pairs).  This is to create a map of the max colors used in
    %plotting (6) and relate it to all the unique sibling pair colors
    color_matrix(1:max(sib_matrix(:,5)),1) = 1:max(sib_matrix(:,5));

    % Creates the map for assigning each sibling color pair a color from 1 to 6


    color_var = 1:6;
    while k <= length(color_matrix) % changed from k = length... 4/17/12
        if color_index == 7
            color_index = 1;
        end
        color_matrix(k,2) = color_var(color_index);
        color_index = color_index+1;
        k = k+1;
    end

    % Map back from the sibling matrix to the output matrix, assigning the
    % unique sibling colors to all cell ids (color assigned after the frame
    % which division begins)
    for s = 1:size(sib_matrix,1);

        %output_index is the mapping index from sibling matrix to xyzs_id
        output_index = sib_matrix(s,1);

        %color_store is the color assigned in the sib_matrix that will be
        %transferred to the xyzs_id matrix
        color_store = sib_matrix(s,5);

        %this loop assigns the color from the sib_matrix to the xyzs_id matrix
        %for all frames after a cell begins to divide
        cell_id = sib_matrix(s,2);
        while cell_id==sib_matrix(s,2)
            xyzs_id(output_index,xyzs_id_columns+1) = color_store;
            output_index = output_index + 1;
            if output_index < length(xyzs_id);
                cell_id = xyzs_id(output_index,xyzs_id_columns);
            else cell_id = 0;
            end;
        end;
    end

    %this loop creates a new vector in the xyzs_id matrix that converts each
    %sibling color to the plotting color (1-to-6) assigned via the color_matrix
    %map.
    xyzs_id(:,xyzs_id_columns+2) = 0;
    for r = 1:length(xyzs_id)
        if xyzs_id(r,xyzs_id_columns+1)>0
            xyzs_id(r,xyzs_id_columns+2) = color_matrix(xyzs_id(r,xyzs_id_columns+1),2);
        end    
    end
    display('Completed sibling sorting')
end

end

