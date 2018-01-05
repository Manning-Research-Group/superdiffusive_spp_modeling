function [] = print_ellipses(pel, Style)
% Function to plot ellipses
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
%  pel: ellipse information for each particle
%  Style: how to plot the ellipse
%
%  OUTPUTS: none
%
%  LOCAL PARAMETERS:

if nargin < 2
    Style = 'r-';
end

for i=1:length(pel)
    hold on;
    es=pel{i};
    if es(1,:) == 0
        disp('particle is a ring')
        disp(i)
    else
    plot(es(1,:), es(2,:), Style);
    end
end

end