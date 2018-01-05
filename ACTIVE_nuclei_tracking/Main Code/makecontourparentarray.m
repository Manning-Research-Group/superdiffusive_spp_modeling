function [ parent_info, parent_vec, x_vector, y_vector ] = makecontourparentarray( nlevels, A, halfobjectsize, noise_wavelength, plottoggle )
% function [ parent_info, parent_vec, x_vector, y_vector, index_values ] = makecontourparentarray( nlevels, A, halfobjectsize, noise_wavelength, plottoggle )
% Function to create a nested list of contours from the image
% intensity_matrix
%
%  1/11/2012
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
%  nlevels: number of contour heights or "levels"
%  A: matrix of intensity values from an image
%  noise_wavelength: characteristic lengthscale of noise in pixels; may be 
%    any positive, floating value.
%  halfobjectsize: A pixel length indicative of a little larger than half a
%    typical object size. 
%  plottoggle – turns plotting capabilities on/off. If set to 1, a contour 
%    plot and a color map (imagesc) are displayed (color map includes peak 
%    contours plotted as white lines). If 0, no plots displayed. 
%
%  OUTPUTS:
%  parent_info: Array defining individual cells. Each entry consists of a
%    set of contours comprising a single cell. 
%  parent_vec: 
%  x_vector: An array of x positions for each contour
%  y_vector: An array of y positions for each contour


% Run bandpass filter on the image
S = bpass(A, noise_wavelength, halfobjectsize); 
intensity_matrix = S;

% Plot the contour map if desired, otherwise just store the contour
% information
if plottoggle ==1;
    figure;
    c = contour(intensity_matrix, nlevels);
else 
    c = contourc(intensity_matrix, nlevels);
end

[mat_row, mat_col] = size(c);
read_index = 1;
num_contours = 0; % actual number of contours lines in c
max_contours = ceil(mat_col/4); %maximum possible number of contours in c

x_vector = cell(max_contours,1); % x_vector{i} is the array of x positions for the ith contour
y_vector = cell(max_contours,1); % y_vector{i} is the array of x positions for the ith contour
c_val = zeros(max_contours,1);  % c_val(i) is the value of the contour
c_numv = zeros(max_contours,1); %c_numv(i) is the number of vertices in each contour
c_COM = zeros(max_contours,2); %center of mass for each contour
c_maxR =zeros(max_contours,1); % maximum radius for each contour

%loop through to save information from contour plot
while (read_index < mat_col) % while the starting point is smaller than the total number of columns
    numv = c(2,read_index);
    value = c(1, read_index);
    
    x_val = c(1, (read_index+1):(read_index+numv));
    y_val = c(2, (read_index+1):(read_index+numv));
    num_contours = num_contours+1;
    
    x_vector{num_contours} = x_val;
    y_vector{num_contours} = y_val;
    
    %calculate the x and y COM for each contour
    c_COM(num_contours,1) =  mean(x_vector{num_contours}); % xCOM
    c_COM(num_contours,2) =  mean(y_vector{num_contours}); % yCOM
    
    %calculate the maximum distance b/w COM and vertices for each contour
    eDist = sqrt((c_COM(num_contours,1)-x_vector{num_contours}).^2+(c_COM(num_contours,2)-y_vector{num_contours}).^2); %calculates Euclidean distance betwen COM and vertices
    c_maxR(num_contours) = max(eDist); %finds maximum
    
    
    c_val(num_contours) = value;
    c_numv(num_contours) = numv;
    %increment index to the next starting point
    read_index = read_index+numv+1;
        
end

% values: vector of unique contour values (repetitions omitted)
% last_pos: vector containing index number of last position for each unique contour value
% index_values: vector of integer level numbers for each contour
[values, last_pos, index_values] = unique(c_val(1:num_contours,:));

%give each contour with level = 1, parent -1
first_pos = find(c_val(1:num_contours,:) > 0, 1); %indentifies first positive contour value

% Valleys represent potential indexing issues - the indexing system is
% re-sorted to accommodate (shouldn't be an issue with cell images)
index_values = index_values-index_values(first_pos)+1;
parent_vec = zeros(num_contours, 3); % (:,1) = index of parent; (:,2) = number of children that parent has, (:,3) = sister flag
index_values = [index_values; NaN]; %break for parent/child while loop

% Added 4/16/2012 - Assign value of -1 for base contours (index values <=1)
bool_par_id = 1;    % lowest level contour
lowest_level = index_values==bool_par_id; %find where lowest level contour occurs
parent_vec(lowest_level,1) = -1; % Assign a value of -1 to lowest level contour (denotes parent level)

j = find(parent_vec(:,1) > -1, 1); %indentifies starting location for j
c_number = [1:num_contours]';   %creates list of contours

%loop to identify parent/child relationships
while index_values(j) <=  max(index_values);
    
    %establish parameters of interest for current level
    COMx = c_COM(j, 1); 
    COMy = c_COM(j, 2);
    start = find(index_values(:,1) == index_values(j)-1,1); % Find first index value of parent level
    ending = find(index_values(:,1) == index_values(j)-1, 1, 'last'); % Find last index value of parent level
    s = 1;
    p_vec_max = length(find(index_values == index_values(j)-1)); % Finds total number of parent vectors
    
    p_COMx = zeros(p_vec_max,1);    % Establish potential parent COM_x vector
    p_COMy = zeros(p_vec_max,1);    % Establish potential parent COM_y vector
    p_maxR = zeros(p_vec_max,1);    % Establish potential parent maxR vector
    p_nums = zeros(p_vec_max,1);    % Find number of potential parents
    
    %find potential parent parameters
    for k = start:ending
        p_COMx(s) = c_COM(k,1);
        p_COMy(s) = c_COM(k,2);
        p_maxR(s) = c_maxR(k);
        p_nums(s) = c_number(k);
        s = s+1;
    end
    
    %identify actual parent using COM and maxR
    eDist = sqrt((COMx-p_COMx).^2+(COMy-p_COMy).^2); %calculates Euclidean distance betwen COM and vertices
    min_Dist = min(eDist);
    myparent = p_nums(eDist == min_Dist);
    parent_vec(j,1) = myparent; %id of the parent
    parent_vec(myparent,2) = parent_vec(myparent,2) +1; % number of children
    
    %ACTION POINT 1/18/2012
    %should we check here to make sure that the contour is actually inside
    %the other (i.e. the max and min x-distances of the interior are less
    %than the max and min x-distances of the parent, and same for y?)
    
    j = j+1;
end


% plot to see if the top of the peaks look correct
if plottoggle == 1;
    figure;
    imagesc(intensity_matrix);
    hold on;
end
peak_index = setdiff(1:num_contours, parent_vec(:,1));
if plottoggle ==1;
    for i=1:length(peak_index)
      plot(x_vector{peak_index(i)},y_vector{peak_index(i)}, 'w-', 'LineWidth', 4);
    end
end


%make a new array to hold nested parent information
% each cell entry contains a vector that starts with a value in peak_index and contains
% all the contours nested below that peak
parent_info = cell(length(peak_index),1); 
for i = 1:length(peak_index)
    current_list = peak_index(i);
    last_parent = parent_vec(peak_index(i),1);

    while last_parent ~= -1
        current_list = [current_list last_parent];
        last_parent = parent_vec(last_parent,1);
    end

    parent_info{i}= current_list;
end

%loop to identify sister contours
k=1;
parent_vec(num_contours+1,:) = NaN; %Added to exit loop
while index_values(k) <=  max(index_values)
    
    current_COMx = c_COM(k, 1);     %x_COM of current contour
    current_COMy = c_COM(k, 2);     %y_COM of current contour
    current_maxR = c_maxR(k);       %maxR of current contour
    sv = find(index_values(:,1) == index_values(k),1);  %starting value - overall index
    % of where this contour level begins
    ev = find(index_values(:,1) == index_values(k), 1, 'last'); %ending value - overall index
    % of where this contour level ends
    
    if ev-sv > 0    %check that more than 1 value exists in this level
    s = 1;  %loop counter
    s_vec_max = length(find(index_values == index_values(k))); %max number of sisters at this level 

    s_COMx = zeros(s_vec_max,1);        %establish sister COM_x vector
    s_COMy = zeros(s_vec_max,1);        %establish sister COM_y vector
    s_maxR = zeros(s_vec_max,1);        %establish sister maxR vector
    s_nums = zeros(s_vec_max,1);        %establish sister number vector
    
    %find potential sister parameters
    for p = sv:ev;
        s_COMx(s) = c_COM(p,1);
        s_COMy(s) = c_COM(p,2);
        s_maxR(s) = c_maxR(p);
        s_nums(s) = c_number(p);
        s = s+1;
    end
    
    %identify sister using COM and maxR
    eDist = sqrt((current_COMx-s_COMx).^2+(current_COMy-s_COMy).^2); %calculates Euclidean distance betwen COM and vertices
    if index_values(k) == 1
        eDist(k) = [];
    else
        num_below = find(index_values(:,1) == index_values(k)-1, 1, 'last');
        eDist(k - num_below) = [];
    end
    
    % establish minimum distance and location
    [minDist minDistloc] = min(eDist);
    
    %Adjust location to accommodate same contour removal (eDist = 0)
    if minDistloc > k
        if minDist < current_maxR && current_maxR > s_maxR(minDistloc + 1)
        parent_vec(minDistloc+1, 3) = 1;
        end
    else
        if minDist < current_maxR && current_maxR > s_maxR(minDistloc)
        parent_vec(minDistloc, 3) = 1;
        end 
    end
    
    end
    k = k+1;
end

%store information in parent_vec
parent_vec(num_contours+1,:) = [];

end
