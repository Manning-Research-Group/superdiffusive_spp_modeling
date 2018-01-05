function [total_cost cost_vector] = cost(xyzs_id1, xyzs_id2, frame_avg, xyzs_id_columns )
% Function to calculate the cost of switching cell IDs during a collision
% event. The cost is determined from cell characteristics that are
% normalized by the frame averages, and then weighting factors allow the
% user to define which parameters are more revealing of a cell's identity.
% 
%  11/10/2012
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
%  xyzs_id1, xyzs_id2: vector of information from xyzs_id for the two
%                      particles to be compared (one from the first frame,
%                      one from the last frame of the collision)
%  frame_avg: matrix of average characteristics of all particles in a given
%             frame (each row is a characteristic and each column is a frame; obtained from frame_char.m)
%
%  OUTPUTS:
%  total_cost: resulting cost from the comparison
%  cost_vector: vector breaking down each component of the final cost
%
%  LOCAL PARAMETERS:
% Determines the weight of each cost parameter (default at 1; can be added
% as user input later)
w_int_intensity = 1;
w_norm_intensity = 1;
w_area = 1;
w_position = 0;
w_aspect = 1;

frame1 = xyzs_id1(xyzs_id_columns-1);
frame2 = xyzs_id2(xyzs_id_columns-1);

% Extract frame averages from frame_avg
avg_norm_intensity1 = frame_avg(7,frame1); %average intensity/area
avg_norm_intensity2 = frame_avg(7,frame2);
avg_area1 = frame_avg(5,frame1); %average area
avg_area2 = frame_avg(5,frame2);
avg_int_intensity1 = frame_avg(6,frame1); %average integrated intensity
avg_int_intensity2 = frame_avg(6,frame2);
avg_diameter1 = (2*frame_avg(1,frame1) + 2*frame_avg(2,frame1))/2; % average of major and minor diameter 
avg_diameter2 = (2*frame_avg(1,frame2) + 2*frame_avg(2,frame2))/2;


% Extracts the intensity
int_intensity1 = xyzs_id1(7);
int_intensity2 = xyzs_id2(7);

% Extracts the areas
area1 = xyzs_id1(6);
area2 = xyzs_id2(6);

% Calculate the average intensity
norm_intensity1 = int_intensity1/area1;
norm_intensity2 = int_intensity2/area2;

% Calculate the aspect ratio
aspect1 = xyzs_id1(3)/xyzs_id1(4);
aspect2 = xyzs_id2(3)/xyzs_id2(4);

% Calculate the distance between centroids
centroid_distance = ((xyzs_id1(1) - xyzs_id2(1))^2 + (xyzs_id1(2) - xyzs_id2(2))^2)^(1/2); 

% Calculate the total cost based on weights and differences
cost_vector(1) = w_int_intensity*(int_intensity1/avg_int_intensity1-int_intensity2/avg_int_intensity2)^2;
cost_vector(2) = w_norm_intensity*(norm_intensity1/avg_norm_intensity1-norm_intensity2/avg_norm_intensity2)^2; 
cost_vector(3) = w_area*(area1/avg_area1 - area2/avg_area2)^2;
cost_vector(4) = w_aspect*(aspect1 - aspect2)^2;
cost_vector(5) = w_position*(centroid_distance/(avg_diameter1/2 + avg_diameter2/2))^2;
total_cost = sum(cost_vector);
cost_vector(6) = total_cost;

end