function [eplotdata] = plot_ellipse(xc, yc, a, b, theta, npoints)
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
%  xc and yc: x and y coordinate, respectively, of ellipse centroid 
%  a and b: major and minor axis, respectively
%  theta: angle of major axis
%  npoints: number of points used to generate plotted ellipse
%
%  OUTPUTS:
%  eplotdata: matrix containing x and y coordinates for an ellipse to be
%  plotted
%  LOCAL PARAMETERS:

% Use equation for an ellipse to determine x and y coordinates that can be
% used to generate a plot of the ellipse
phi = 2*pi*(1:npoints)/npoints;
eplotdata = zeros(2, npoints);
eplotdata(1,:) = xc + a*cos(theta)*cos(phi) - b*sin(theta)*sin(phi);
eplotdata(2,:) = yc + a*sin(theta)*cos(phi) + b*cos(theta)*sin(phi);

end
