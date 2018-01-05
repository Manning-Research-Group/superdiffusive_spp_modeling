function [frame_avg] = frame_char(xyzs_id,xyzs_id_columns)
% function [frame_charact, frame_avg, frame_std] = frame_char(xyzs_id,xyzs_id_columns)
%
% Function to determine frame averages of the major axis, minor axis,
% aspect ratio, orientation angle, area, integrated intensity and average
% intensity.
% 
%  10/1/2012
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
%  xyzs_id: matrix of particle information of tracked particles (from
%           trackmem_new.m)
%  xyzs_id_columns: column in xyzs_id that contains cell ID information
%
%  OUTPUTS:
%  frame_avg: matrix containing the average of each characteristic for each frame
%  
%  ADDITIONAL OPTIONAL OUTPUTS:
%  frame_charact: cell array containing particle information for each frame
%  frame_std: matrix containing the standard deviation of each characteristic for each frame  
%
%  LOCAL PARAMETERS:

for z = 1:max(xyzs_id(:,xyzs_id_columns-1)); % z corresponds to each frame
    bool_vector = xyzs_id(:,xyzs_id_columns-1) == z; % bool vector for extracting info for each frame
    frame_charact{1,z} = real(xyzs_id(bool_vector,3)); % major axis 
    frame_charact{2,z} = real(xyzs_id(bool_vector,4)); % minor axis 
    frame_charact{3,z} = real(frame_charact{1,z}./frame_charact{2,z}); % aspect ratio 
    frame_charact{4,z} = real(xyzs_id(bool_vector,5)); % angle of orientation 
    frame_charact{5,z} = real(xyzs_id(bool_vector,6)); % area
    frame_charact{6,z} = real(xyzs_id(bool_vector,7)); % integrated intensity
    frame_charact{7,z} = real(frame_charact{6,z}./frame_charact{5,z}); % average intensity
    
    % Calculate averages
    frame_avg(1,z) = mean(frame_charact{1,z});
    frame_avg(2,z) = mean(frame_charact{2,z});
    frame_avg(3,z) = mean(frame_charact{3,z});
    frame_avg(4,z) = mean(frame_charact{4,z});
    frame_avg(5,z) = mean(frame_charact{5,z});
    frame_avg(6,z) = mean(frame_charact{6,z});
    frame_avg(7,z) = mean(frame_charact{7,z});
    
    % Calculate standard deviations
    frame_std(1,z) = std(frame_charact{1,z});
    frame_std(2,z) = std(frame_charact{2,z});
    frame_std(3,z) = std(frame_charact{3,z});
    frame_std(4,z) = std(frame_charact{4,z});
    frame_std(5,z) = std(frame_charact{5,z});
    frame_std(6,z) = std(frame_charact{6,z});
    frame_std(7,z) = std(frame_charact{7,z});
end