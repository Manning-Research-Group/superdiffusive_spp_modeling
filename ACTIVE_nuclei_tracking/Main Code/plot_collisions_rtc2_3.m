function [  ] = plot_collisions_rtc2_3( xyzs_id, event_array, stackname, xyzs_id_columns, data_folder, post_correct, video_disp_toggle)
%   Function that plots videos for division and merging events.
%   For each event detailed in the event array, this function plots a video
%   of the corresponding trajectories for all cells involved in the
%   division or merging event. This allows the user to manually check
%   trajectories to ensure proper cell identification. 
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
%  xyzs_id: matrix of particle information post-tracking and/or
%           post-correcting
%  event_array: cell_array output from collision_corrector.m that contains
%               information for every event
%  stackname: name of the stack being analyzed; used for creating the video
%
%  OUTPUTS:
%  
%  LOCAL PARAMETERS:

if nargin < 6
    post_correct = 0;
end

if nargin < 7
    video_disp_toggle = 1;
end

image_info = imfinfo(stackname);
number_images = numel(image_info);

% Create new directory to store video information
if post_correct == 1;
   videofolder = 'events_post_correcting';
else
    videofolder = 'events_pre_correcting';
end

mkdir(data_folder, videofolder);
type_chk = strfind(stackname, '/');
if ~isempty(type_chk)
    addpath(strcat(data_folder, '\', videofolder));
else
    addpath(strcat(data_folder, '/', videofolder));
end
nume = size(event_array,1); % Number of events
figbase = 'event'; % Figure name base

color_mat = ['b', 'g', 'r', 'y', 'm', 'c']; % Track trajectory colors for video plotting

%In general, we want a subset of the events printed (10-15 events max); 
%when the total event number is smaller, you need to multiply by a larger
%factor to achieve the same number of events. 
if nume>100
    mod_val = ceil(nume*.10);
else
    mod_val = ceil(nume*.30);
end
es = mod([1:nume],mod_val);   %Event spacing
bool_es = es == 1;      %Identify 
output_events = find(bool_es);   %Specify which event videos to output

for i = output_events
    % Each merging event will contain different numbers interacting cells.
    % The cleared variables must be refreshed for each event to ensure
    % proper video window adjustment and cell information. 
    clear xmin xmax ymin ymax max_frame_vec min_frame_vec cell_array 
    clear min_xpos max_xpos min_ypos max_ypos
    event_info = event_array{i};    %Event information for 1 event

    % Identify interacting cells
    unique_cells = nonzeros(unique(event_info(:,2:3)));
    num_cells = length(unique_cells);
    cell_array = cell(num_cells,1);

    % For each cell, find its merging event xyzs_id information and identify
    % a video window for viewing
    for k = 1:num_cells
        clear cell_mat sib_mat cell_sort sib_sort frame
        cell_pos = xyzs_id(:,xyzs_id_columns) == unique_cells(k);
        cell_mat = xyzs_id(cell_pos,:);
        cell_sort = sortrows(cell_mat,xyzs_id_columns-1);
        cell_array{k} = cell_sort;
        
        xpos{k} = cell_sort(:,1);
        min_xpos(k) = floor(min(cell_sort(:,1)));
        max_xpos(k) = ceil(max(cell_sort(:,1)));
        ypos{k} = cell_sort(:,2);
        min_ypos(k) = floor(min(cell_sort(:,2)));
        max_ypos(k) = ceil(max(cell_sort(:,2)));
        
        max_frame_vec(k) = max(cell_sort(:,xyzs_id_columns-1));
        min_frame_vec(k) = min(cell_sort(:,xyzs_id_columns-1));
        
    end

    % Identify overall video window for all interacting cells
    xmin = min(min_xpos)-30;
    ymin = min(min_ypos)-30;
    xmax = max(max_xpos)+30;
    ymax = max(max_ypos)+30;
    
    % If clauses to prevent extraneous window formatting 
    if xmin < 1;
        xmin = 1;
    end
    
    if xmax > image_info(1).Width;
        xmax = image_info(1).Width;
    end

    if ymin < 1;
        ymin = 1;
    end
    
    if ymax > image_info(1).Height;
        ymax = image_info(1).Height;
    end
    
    % Identify overall frame information
    max_frame = max(max_frame_vec);
    min_frame = min(min_frame_vec);

    % Find first frame of interaction
    current_frame = min_frame;
 
    %Determine figure name (Event#)
    figname = [figbase, '_',num2str(i)];
    
    if video_disp_toggle == 1
        % Set-up video writer
        frame(1:size(cell_mat,1)) = struct('cdata', [],...
                'colormap', []);
            
%         cd([save_path '\' videofolder])
        if ~isempty(type_chk)
            writerObj = VideoWriter([data_folder '\' videofolder '\' figname]);
        else
            writerObj = VideoWriter([data_folder '/' videofolder '/' figname]);
        end
        writerObj.FrameRate = 5;
        open(writerObj);
%         cd ..
    else
        %Make directory (event#)
%         cd([save_path '\' videofolder])
        
        if ~isempty(type_chk)
            mkdir([data_folder '\' videofolder], figname);
            addpath([data_folder '\' videofolder '\' figname]);
        else
            mkdir([data_folder '/' videofolder], figname);
            addpath([data_folder '/' videofolder '/' figname]);
        end
%         cd ..
    end
    
    % Loop through all frames and plot corresponding cell merging
    % information. 
    if max_frame > min_frame
    for r = 1:(max_frame - min_frame);
        
        % Start by plotting full image and incorporating cell ID values
        h = figure(i);
        hold off;
        A = imread(stackname, current_frame, 'Info', image_info);
        imagesc(A);
        colormap('gray');
        hold on
        bool_frame = xyzs_id(:,xyzs_id_columns-1)==current_frame;
        text(xyzs_id(bool_frame,1),xyzs_id(bool_frame,2), num2str(xyzs_id(bool_frame,xyzs_id_columns)),'FontSize', 14, 'Color', 'y');
        
        % Change to merging event window
        axis([xmin xmax ymin ymax])
        
        % Create array for holding cell index information (from xyzs_id)
        clear bool
        bool = cell(1,num_cells);
        
        % Loop to find and plot information for all cells involved in
        % merging event
        for j = 1:num_cells
            bool{j} = cell_array{j}(:,xyzs_id_columns-1)<=current_frame;
            plot(cell_array{j}(bool{j},1), cell_array{j}(bool{j},2),'marker', 'o', 'color', color_mat(j));
            hold on
        end

        % Move to next frame and repeat process
        current_frame = current_frame + 1;
        
        if video_disp_toggle == 1
            % Advance video editor frame
            frame(r) = getframe;
        else
            real_figname = ['event', num2str(i), '_', num2str(r)];
            if ~isempty(type_chk)
                saveas(h, [data_folder '\' videofolder '\' figname '\' real_figname] , 'tif');
            else
                saveas(h, [data_folder '/' videofolder '/' figname '/' real_figname] , 'tif');
            end
            close(h);
        end
    end
    
        if video_disp_toggle == 1
            % Save and write information to video file
            writeVideo(writerObj,frame(1:r));
            close(writerObj);
            close(h);
        end
    end

fprintf('\nEvent %d complete\n',i)

end