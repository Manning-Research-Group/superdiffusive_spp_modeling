function [] = analysis2(angle_vec, plot_toggle)
%The analysis function is used to generate correlation plots for multiple
%samples. 
    % Inputs: 
    %   1) avg_tog - toggle to caculate sample averages. 1 = on, 0 = off
    %   2) angle_vec - row vector of angle values (radians) for anisotropic 
    %   samples (used for axis shiftin); this vector will default to 0's if 
    %   no values are specified. Angle gvalues should fall between -2pi and
    %   2pi
    %   2) accuracy - accepted values include 2, 4, 6, or 8. Defines the
    %   number of points used in the central finite difference calculation
    %   for velocity
    %   3) framerate - time lapse (minutes) between two consecutive frames
    %
    % Outputs: 
    %   MSD Plots: 1) SampleName_full_MSD: MSD plot for overall mean
    %   squared displacement calculation with polynomial fit and
    %   upper/lower slopes descriptor; 2) SampleName_full_alpha: plot of
    %   non-Gaussian parameter for overall MSD calculation
    %   
    %   MSD Decomposition Plots: 1) SampleName_tracks: plots trajectories
    %   for all tracked cells over time; 2) SampleName_alpha_no_fit: plots
    %   overall, x, and y non-Gaussian parameters in a single plot; 3)
    %   SampleName_alpha_power_fit: same as SampleName_alpha_no_fit but
    %   with a power law fit to each curve; 4) SampleName_MSD_decomp: MSD
    %   components for x and y directions with polynomial fit and
    %   upper/lower slopes descriptors. 
    %   
    %   Velocity Autocorrelation Plots: 1) SampleName_velocity: plot of x
    %   and y velocities with exponential fit. 
    %
    %   .mat files: 1) MSD_decomp: x,y data pertaining to decomposed 
    %   non-Gaussian parameter; 2) velocity.mat: updated xyzs_id matrix 
    %   with velocity calculations and x,y information from velocity plot. 

    % Megan Brasch, Richard Baker 6/4/13
accuracy = 8;
framerate = 3;
    
if nargin < 2;
    plot_toggle = 1;
end

%Find all of the .mat files present in the current directory
file_info = dir('*.mat');
m = size(file_info,1);

%Make angle_vec if one is not input
if nargin<1
    angle_vec = zeros(m,1);
end

%If partial angle_vec is input, add zeros for all other recognized samples
if size(angle_vec,2)~=m
    diff = m-size(angle_vec,2);
    extra = zeros(1,diff);
    angle_vec = [angle_vec, extra]; %Non-wrinkled or TCPS samples only (alphabetical)
end

original_dir = pwd;

%Find all angles that are not zero
angle_vec_idx = angle_vec ~= 0;
angle_vec_chk = angle_vec(angle_vec_idx);

orig_angles = angle_vec;
%Convert angles to radians if input as degrees
if (angle_vec_chk(:) > 2*pi) | (angle_vec_chk(:) < -2*pi)
    angle_vec = angle_vec.*(pi/180);
end

%Make matrices to store MSD, decomposed MSD, and velocity constants
MSD_mat = zeros(m,3);   %columns: 1) short timescale slope, 2) long timescale slope, 3) difference of short and long timescale slopes
decomp_MSD_mat = zeros(m,7); 
vel_mat = zeros(m,6);

%Arrange output list for Excel file writing
output_list = cell(m, 33);
output_list{1,1} = 'Filename';
%Full MSD ouputs:
output_list{1,2} = 'Short Timescale Slope';
output_list{1,3} = 'Long Timescale Slope';
output_list{1,4} = 'Long Timescale Intercept';
output_list{1,5} = 'Difference (short-long)';
output_list{1,6} = 'Intersect X Coord';
output_list{1,7} = 'Intersect Y Coord';
%Decomposed MSD outputs:
output_list{1,8} = 'X Short Timescale Slope';
output_list{1,9} = 'X Long Timescale Slope';
output_list{1,10} = 'X Long Timescale Intercept';
output_list{1,11} = 'Y Short Timescale Slope';
output_list{1,12} = 'Y Long Timescale Slope';
output_list{1,13} = 'Y Long Timescale Intercept';
output_list{1,14} = 'X Difference (short-long)';
output_list{1,15} = 'Y Difference (short-long)';
output_list{1,16} = 'Difference of differences (x diff-y diff)';
output_list{1,17} = 'x_Intersect X Coord';
output_list{1,18} = 'x_Intersect Y Coord';
output_list{1,19} = 'y_Intersect X Coord';
output_list{1,20} = 'y_Intersect Y Coord';
output_list{1,21} = 'Diff X Coord (x-y)';
output_list{1,22} = 'Diff Y Coord (x-y)';
%Velocity outputs: 
output_list{1,23} = 'X Fit A Constant';
output_list{1,24} = 'X Time Constant (B)';
output_list{1,25} = 'X Fit C Constant';
output_list{1,26} = 'Y Fit A Constant';
output_list{1,27} = 'Y Time Constant (B)';
output_list{1,28} = 'Y Fit C Constant';
output_list{1,29} = 'Time Constant Difference (x-y)';
%Cell Count outputs
output_list{1,30} = 'Number Cells First Frame';
output_list{1,31} = 'Number Cells Last Frame';
output_list{1,32} = 'Number Cells Average';
output_list{1,33} = 'Total Number Cells';

%Make a new directory for each sample (.mat file)
for i = 1:m
    %Identify sample and make new directory to store analysis information
    full_filename = file_info(i).name;
    [~, filename] = fileparts(full_filename);  %Need filename without extension to make directory
    output_list{i+1,1} = filename;
    new_dir = [original_dir, '\', filename];
    load(full_filename)
    mkdir(new_dir)
    theta = angle_vec(i);
    
    %Full MSD analysis
    [short_slope long_slope, intersect, intercept_upper] = MSD_henderson(xyzs_id, new_dir, filename, plot_toggle); 
    output_list{i+1,2} = short_slope;
    output_list{i+1,3} = long_slope;
    output_list{i+1,4} = intercept_upper;
    output_list{i+1,5} = short_slope-long_slope;
    output_list{i+1,6} = intersect(1);
    output_list{i+1,7} = intersect(2);
    MSD_mat(i,1) = short_slope; MSD_mat(i,2) = long_slope; MSD_mat(i,3) = short_slope-long_slope;
    
    %MSD Decomposition
    [xshort_slope, xlong_slope, yshort_slope, ylong_slope, alpha2x, alpha2y, alpha2b, tMSD, intersect_x, intersect_y, intercept_upper_x, intercept_upper_y] = MSD_henderson_theta(xyzs_id, theta, new_dir, filename, plot_toggle);
    save([new_dir, '\MSD_decomp.mat'], 'alpha2x', 'alpha2y', 'alpha2b', 'tMSD', 'intersect_x', 'intersect_y','intercept_upper_x','intercept_upper_y');
    output_list{i+1,8} = xshort_slope;
    output_list{i+1,9} = xlong_slope;
    output_list{i+1,10} = intercept_upper_x;
    output_list{i+1,11} = yshort_slope;
    output_list{i+1,12} = ylong_slope;
    output_list{i+1,13} = intercept_upper_y;
    output_list{i+1,14} = xshort_slope-xlong_slope;
    output_list{i+1,15} = yshort_slope-ylong_slope;
    output_list{i+1,16} = output_list{i+1,9}-output_list{i+1,10};
    output_list{i+1,17} = intersect_x(1);
    output_list{i+1,18} = intersect_x(2);
    output_list{i+1,19} = intersect_y(1);
    output_list{i+1,20} = intersect_y(2);
    output_list{i+1,21} = intersect_x(1)-intersect_y(1);
    output_list{i+1,22} = intersect_x(2)-intersect_y(2);
    decomp_MSD_mat(i,1) = xshort_slope; decomp_MSD_mat(i,2) = xlong_slope; decomp_MSD_mat(i,3) = yshort_slope; decomp_MSD_mat(i,4) = ylong_slope; decomp_MSD_mat(i,5) = xshort_slope-xlong_slope; decomp_MSD_mat(i,6)= yshort_slope-ylong_slope; decomp_MSD_mat(i,7) = decomp_MSD_mat(i,5)-decomp_MSD_mat(i,6);
    
%     %Velocity Autocorrelation
%     [xyzs_id] = velocity_calc2(xyzs_id, xyzs_id_columns, accuracy);
%     [ xa xtime xc ya ytime yc Mcx Mcy NMSD timewindows fstart xv yv] = velocity_auto_theta( xyzs_id, xyzs_id_columns, new_dir, filename, theta, framerate, plot_toggle);
%     save([new_dir, '\velocity.mat'], 'xyzs_id', 'Mcx', 'Mcy', 'NMSD', 'timewindows');
%     output_list{i+1,23} = xa;
%     output_list{i+1,24} = xtime;
%     output_list{i+1,25} = xc;
%     output_list{i+1,26} = ya;
%     output_list{i+1,27} = ytime;
%     output_list{i+1,28} = yc;
%     output_list{i+1,29} = xtime-ytime;
%     vel_mat(i,1) = xa; vel_mat(i,2) = xtime; vel_mat(i,3) = xc; vel_mat(i,4) = ya; vel_mat(i,5) = ytime; vel_mat(i,6) = yc; vel_mat(i,7) = xtime-ytime;
%     
%     %Cell Number information
%     [cell_count_first, cell_count_last, total_cells] = num_cells(xyzs_id, xyzs_id_columns, n_divisions);
%     output_list{i+1,30} = cell_count_first;
%     output_list{i+1,31} = cell_count_last;
%     output_list{i+1,32} = (cell_count_first+cell_count_last)/2;
%     output_list{i+1,33} = total_cells;
%     
   % Diffusion plotting
    if plot_toggle == 1;
        diffusion_plot2(xyzs_id,xyzs_id_columns, theta,new_dir,filename);
    end
%     
% [lambda1, lambda2, radius_gyration, track_length, gyration_angle{i}] = gyration(xyzs_id,xyzs_id_columns);
%     if plot_toggle == 1;
%         gyration_plot2(xyzs_id,xyzs_id_columns, gyration_angle{i},new_dir,filename);    
%     end
    
    close all
end

save([original_dir, '\','Final_Analysis.mat'], 'output_list', 'MSD_mat', 'decomp_MSD_mat', 'vel_mat', 'orig_angles', 'angle_vec');

%Write output to csv file - hopefully Mac compatible!
fid=fopen([original_dir, '\', 'Analysis_Summary.csv'],'wt');

[rows,cols]=size(output_list);

for j=1:rows
    if j == 1
        fprintf(fid,'%s,', output_list{j, 1:end-1});
        fprintf(fid,'%s\n',output_list{j,end});
    else
        fprintf(fid,'%s,', output_list{j, 1});
        fprintf(fid,'%3.4d,', output_list{j,2:end-1});
        fprintf(fid,'%3.4d\n',output_list{j,end});
    end
end

fclose all;

end