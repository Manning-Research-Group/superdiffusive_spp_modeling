function [ xyzs_id, xyzs_id_columns, filename, framerate, new_dir] = run_tracking_contour2(image_mat_in, movie_filename,inputfilename)
% Function to find features and track them from a tiff stack of input data
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
%  inputfilename: location of tiff stack image
%
%
%  OUTPUTS:
%  
%  xyzs_id: 1) x-position, 2) y-position, 3) major axis, 4) minor axis, 5)
%  angle (theta), 6) area, 7) intensity, 8) sibling tag (within frame), 9)
%  sibling tag (overall frames), 10) overall individual cell tag (for 
%  sibling identification), 11) frame number, 12) cell ID, 13) color
%  identification (for movie generation), 14) color identifier scaled to
%  pre-set Matlab colors (for movie generation). 
%
%  sib_matrix: PRE-colision correction: 1) cell 1 index into xyzs_id, 2)
%  cell 1 ID, 3) cell 2 ID, 4) frame number, 5) color number reference
%  (relates to xyzs_id column 13 and 14). POST-collision correction: 6)
%  Type of collision (division, collision, continuing collision), 7) index
%  for event array building, 8) event number, 9) cell 2 index to xyzs_id


%  LOCAL PARAMETERS:
frust_toggle = 1; % Can change this for future experiments without frustrated cell divisions
fit_height = 1/2; % contour level to fit an ellipse to
removeflagged=1; % removes particles that are too small when set to 1

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in user input file if one is supplied, otherwise prompt user for
% parameters
if nargin < 3
    prompt = {'Plot Toggle (1=on; 0=off):','Number of Contours:','Half Particle Diameter (must be an odd integer)','Noise Wavelength (pixels):','Collision Plot Toggle (1=yes; 0=no)','Maximum Area (pixels^2)', 'Minimum Area (pixels^2)', 'Maximum Displacement (pixels)','Frame Time (min)','Maximum Collision Time (frames)', 'Manual Division GUI Toggle (1=yes; 0=no)', 'Manual Merging GUI Toggle (1=yes; 0=no)'};
    dlg_title = 'Tracking Parameters';
    num_lines = 1;
    def = {'0','10','13','2','0','260', '10', '20','3','10', '0', '0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);

    plottoggle = str2num(answer{1});
    nlevels = str2num(answer{2});
    halfobjectsize = str2num(answer{3});
    noise_wavelength = str2num(answer{4});
    collision_plot_toggle = str2num(answer{5});
    area_thresh = str2num(answer{6});
    min_area = str2num(answer{7});
    maxdisp = str2num(answer{8});
    framerate = str2num(answer{9});
    max_collision_time = str2num(answer{10});
    div_toggle = str2num(answer{11});
    merging_toggle = str2num(answer{12});
       
    
else
    fid = fopen(inputfilename);
    values = textscan(fid, '%f %*[^\n]');
    name = textscan(fid, '%s %*s');
    full_name = name{1}{1};
    
    nlevels = values{1}(1);
    halfobjectsize = values{1}(2);
    noise_wavelength = values{1}(3);
    collision_plot_toggle = values{1}(4);
    area_thresh = values{1}(5);
    min_area = values{1}(6);
    maxdisp = values{1}(7);
    framerate = values{1}(8);
    max_collision_time = values{1}(9);
    
    %Command line does not support graphics interface - default any
    %plotting parameters to 0
    plottoggle = 0;
    video_disp_toggle = 0;
    merging_toggle = 0;
    div_toggle = 0;
end;

%allow user to select tif image or tif stack to be loaded
if nargin < 2
    [filename, pathname] = uigetfile({'*.tif'}, 'Select a tif file');
    full_name = [pathname, filename];
    
    video_disp_toggle = 1;
else 
    full_name = movie_filename;
end

[pathname, filename, ext] = fileparts(full_name);
close all
type_chk = strfind(pathname, '\');
if ~isempty(type_chk)
    idx = max(strfind(pathname,'\'));
    home_dir = pathname(1:idx);
    output_path = [home_dir, 'Analysis'];
    mkdir(output_path);
    data_folder = [output_path, '\', filename];
    mkdir(data_folder);
    addpath(pathname,home_dir, output_path, data_folder);
    stackname = [pathname, '\', filename, ext];
else
    idx = max(strfind(pathname,'/'));
    home_dir = pathname(1:idx);
    output_path = [home_dir, 'Analysis'];
    mkdir(output_path);
    data_folder = [output_path, '/', filename];
    mkdir(data_folder);
    addpath(pathname,home_dir, output_path, data_folder);
    stackname = [pathname, '/', filename, ext];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

image_info = imfinfo(stackname);
number_images = numel(image_info);
weird_vec_cell = cell(number_images,1);
fprintf('\nStarting contour analysis and ellipse fitting for each frame\n');
toc;

particle_array = cell(number_images,1); 
pellipses_array = cell(number_images,1);
for j = 1:number_images;
    %send each one of the images in tiff stack to particle finding algorithm
    %save particle information in large cell array

    A = imread(stackname, j, 'Info', image_info); % Read in image data
    if size(image_mat_in,1) == 1
        A = imadjust(A, image_mat_in,[0;1]);
    end
    if size(image_mat_in,1) == number_images
        A = imadjust(A, image_mat_in(j,:),[0;1]);
    end
    [ parent_info, parent_vec, x_vector, y_vector] = makecontourparentarray( nlevels, A, halfobjectsize, noise_wavelength, plottoggle ); % Run contour analysis on an image
    fprintf('\nContours completed for frame %d \n',j);
    toc;
    
    [particles pellipses weird_vec] = find_particles_fixed( parent_info, parent_vec, x_vector, y_vector, removeflagged, area_thresh, min_area, fit_height ); % Detect particles from contour profile and fit with an ellipse
    fprintf('Ellipses completed for frame %d \n',j);
    toc;
    
    particles = ellipse_mask3(particles, A); % Run masking analysis to calculate sum of intensity and area for each particle
    fprintf('Intensity calculations completed for frame %d \n',j);
    toc;
    
    num_particles = size(particles,1);
    fprintf('Number of particles is %d \n', num_particles);

    weird_vec_cell{j} = weird_vec;
    particle_array{j} = particles; 
    pellipses_array{j} = pellipses;
    
    % Will print ellipse fit if plottoggle is 1
    if plottoggle == 1;
        print_ellipses(pellipses, 'r-');
    end;
    
end;

save('variables.mat','particle_array'); 

%loop through cell array and save particle information in tracking format
t_index = 1;
disp('Creating matrix of particle positions and particle descriptors');
%collect features of interest in large t_matrix to send through tracking
for k = 1:size(particle_array,1);
    t_matrix(t_index:t_index + size(particle_array{k},1) - 1,1) = particle_array{k}(:,1); %x position
    t_matrix(t_index:t_index + size(particle_array{k},1) - 1,2) = particle_array{k}(:,2); %y position
    t_matrix(t_index:t_index + size(particle_array{k},1) - 1,3) = particle_array{k}(:,3); %major axis
    t_matrix(t_index:t_index + size(particle_array{k},1) - 1,4) = particle_array{k}(:,4); %minor axis
    t_matrix(t_index:t_index + size(particle_array{k},1) - 1,5) = particle_array{k}(:,5); %angle info
    t_matrix(t_index:t_index + size(particle_array{k},1) - 1,6) = particle_array{k}(:,size(particles,2)-1); %area info
    t_matrix(t_index:t_index + size(particle_array{k},1) - 1,7) = particle_array{k}(:,size(particles,2)); %integrated intensity info
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%$$ START CHANGES

    t_matrix(t_index:t_index + size(particle_array{k},1) - 1,8) = particle_array{k}(:,10); % multi-body interaction flag

    
    t_matrix(t_index:t_index + size(particle_array{k},1) - 1,9) = particle_array{k}(:,9); %sibling info
    t_matrix(t_index:t_index + size(particle_array{k},1) - 1,12) = k; %frame number
    t_index = t_index + size(particle_array{k},1);
end;
t_matrix(:,11) = (1:size(t_matrix,1)); %Used for referencing between t_matrix
disp('Finished matrix of particle positions and particle descriptors');
toc;

save('variables.mat','t_matrix'); 

frame_index = 1; %initialize frame number
total_particles = 0; %initialize count for particles
p_index = 1; %initialize particle index
previous_frame_total = 0;

% convert frame by frame sibling information to overall cell ID numbers
while frame_index <= size(particle_array,1)
    num_particles = size(particle_array{frame_index},1); %identify total particles in that frame
    total_particles = total_particles + num_particles; %identify running total of particles in all previous frames
    
    if frame_index > 1
        previous_frame_total = total_particles - num_particles;
        t_matrix(p_index:total_particles, 10) = t_matrix(p_index:total_particles, 9) + previous_frame_total;
    else
    t_matrix(p_index:total_particles, 10) = t_matrix(p_index:total_particles, 9);
    end
    
    for i = p_index:total_particles
        if t_matrix(i, 10) == previous_frame_total && t_matrix(i, 9) == 0
            t_matrix(i, 10) = 0;
        end
    end
    %%$$ END CHANGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    p_index = total_particles + 1;
    frame_index = frame_index+1;
end

save('variables.mat','t_matrix'); 

disp('Tracking particles')
toc;
%send tracking formatted data to trackmem
dim = 2;
goodenough = 1;
memory = 10;
[xyzs_id] = trackmem_new(t_matrix,maxdisp,dim,goodenough,memory); % Run tracking analysis on t_matrix
disp('Tracking particles completed');
toc;

xyzs_id_columns = size(xyzs_id,2);
save('variables.mat', 't_matrix', 'xyzs_id', 'xyzs_id_columns'); 

% Extract sibling information
disp('Extracting sibling information');
toc;
[sib_matrix, xyzs_id] = sib_matrix_creation(xyzs_id, xyzs_id_columns);
save('variables.mat', 'xyzs_id', 'sib_matrix', 'xyzs_id_columns','t_matrix');

if ~isempty(sib_matrix)
    % Determine particle characteristics for each frame
    [frame_avg] = frame_char(xyzs_id, xyzs_id_columns);

    [xyzs_id, sib_matrix, event_array xyzs_id_columns] = collision_corrector_rtc2(xyzs_id, sib_matrix, xyzs_id_columns,max_collision_time);

    fprintf('Total number of events: %d\n', size(event_array),1);

    if collision_plot_toggle == 1
        plot_collisions_rtc2_3( xyzs_id, event_array, stackname, xyzs_id_columns, data_folder, 0, video_disp_toggle);
    end

    [xyzs_id, mastertable] = collision_tags_rtc2_3(xyzs_id, event_array, xyzs_id_columns, frame_avg);
    [mult_sib_array, chk_mult_sib] = mult_sib_creation(xyzs_id,xyzs_id_columns);
    if chk_mult_sib == 1
        [mult_array] = mult_event_array(mult_sib_array,max_collision_time);
    else
        mult_array = [];
    end
    
    %Call manual merging GUI based on user input
    if merging_toggle == 1
        [merging_output, pause_output] = merging_gui(stackname, mult_array, xyzs_id, xyzs_id_columns);
    else
        merging_output = [];
        pause_status = [];
    end

    %Send info to new function to combine mastertable and merging_gui output
    [compiled_mastertable] = master_compile(mastertable, merging_output);

    disp('Completed identification of complex merging events.')

    disp('Processing division information...')
    [xyzs_id event_array missing_cell_mat] = mult_division_adjuster(xyzs_id,xyzs_id_columns,event_array, mult_array, max_collision_time);

    [xyzs_id event_array2 div_angle n_divisions] = division_corrector(xyzs_id, event_array, xyzs_id_columns, halfobjectsize, A, div_toggle, 5, 30, frust_toggle, stackname);
    disp('Completed division processing')

    % Perform relabelling analysis
    [ xyzs_id ] = relabel(compiled_mastertable, xyzs_id, xyzs_id_columns);

    xyzs_id = sortrows(xyzs_id, [xyzs_id_columns xyzs_id_columns-1]);

    if collision_plot_toggle == 1
        plot_collisions_rtc2_3( xyzs_id, event_array, stackname, xyzs_id_columns, data_folder, 1, video_disp_toggle);
    end
    
    if ~isempty(type_chk)
        save_name = strcat(data_folder, '\', filename, '.mat');
    else
        save_name = strcat(data_folder, '/', filename, '.mat');
    end

    new_dir = data_folder;
    save(save_name,'xyzs_id', 'sib_matrix', 'event_array', 'xyzs_id_columns', 'frame_avg', 't_matrix', 'event_array2', 'div_angle', 'n_divisions', 'mult_sib_array', 'mult_array', 'new_dir','filename');

else % If no cell interactions just save non-event information
    if ~isempty(type_chk)
        save_name = strcat(data_folder, '\', filename, '.mat');
    else
        save_name = strcat(data_folder, '/', filename, '.mat');
    end

    new_dir = data_folder;
    save(save_name,'xyzs_id', 'xyzs_id_columns', 't_matrix', 'new_dir','filename');
end

% Uncomment the following to create a movie with overlay ellipses
% [a ] = make_ellipse_movie(xyzs_id, xyzs_id_columns, stackname, data_folder);

end