function [particles] = ellipse_mask3(particles, A)
% Function to create an elliptical mask for each particle in an image that
% is used to determine the sum of intensity for that particle; this sum of
% intensity along with the calculated area is added to the particles matrix
% and returned
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
%  particles: matrix of particle information where first two columns must
%             be the x and y coordinates, respectively, of the center of the ellipse
%             and the third and fourth columns must be the major and minor axis,
%             respectively
%  A: matrix of intensity values from an image
%
%  OUTPUTS:
%  particles: same matrix as input with two added columns: 1)area of each
%             particle and 2)sum of intensity for each particle
%
%  LOCAL PARAMETERS:

% Add two columns of zeros to the particles matrix for space allocation
num_columns = size(particles,2);
particles(:,num_columns+1:num_columns+2) = 0;

particles(particles(:,1)<0,:) = [];
particles(particles(:,2)<0,:) = [];
particles(particles(:,1)>size(A,2),:) = [];
particles(particles(:,2)>size(A,1),:) = [];

% Loop to calculate the sum of intensity and area of a particle
for i = 1:size(particles,1)
    xc = particles(i,1);
    yc = particles(i,2);

    a = particles(i,3);
    b = particles(i,4);
    theta = particles(i,5);
    
    % determines the width and height of the smallest matrix that can
    % encapsulate the particle of interest (used to increase masking speed)
    if a*cos(theta) < b
        width = ceil(b);
    else
        width = ceil(a*cos(theta));
    end
    if a*sin(theta) < b
        height = ceil(b);
    else
        height = ceil(a*sin(theta));
    end
  
    % Creates a cropped matrix of the original image data that contains the
    A_prime = A(max(0,ceil(yc)-height):min(size(A,1),ceil(yc)+height), max(0,ceil(xc)-width):min(size(A,2),ceil(xc)+width));
    
    % Determines tranformation matrices for the "center form" equation of an ellipse 
    R_matrix = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    D_matrix = [1/(a^2) 0; 0 1/(b^2)];
    % Center form equation for an ellipse
    A_matrix = R_matrix*D_matrix*R_matrix';

    % Creates the masking matrix where anything within the ellipse is a "1"
    % and everything outside the ellipse is a "0"
    yspan = abs(min(size(A,1),ceil(yc)+height)- max(0,ceil(yc)-height));
    xspan = abs(min(size(A,2),ceil(xc)+width) - max(0,ceil(xc)-width));
    [X,Y] = meshgrid(0:xspan,0:yspan);
    
    X = X - width;
    Y = Y - height;
    mask = X.^2 * A_matrix(1,1) + X.*Y * (A_matrix(1,2) + A_matrix(2,1)) + Y.^2 * A_matrix(2,2) < 1;
    
    % Calculates the area; also calculates the sum of intensity by
    % overlaying the mask to the cropped image and summing
    particles(i,num_columns+1) = pi*particles(i,3)*particles(i,4); % calculate the area for particle
    particles(i,num_columns+2) = sum(A_prime(mask)); % calculate sum of intensity for particle
end