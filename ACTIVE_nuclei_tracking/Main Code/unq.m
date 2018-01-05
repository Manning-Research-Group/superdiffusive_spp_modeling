% Used for Kilfoil linking

%THIS CODE IS PLANNED FOR ARRAY BEING A ROW VECTOR
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

function [ret] = unq(array,idx)

s = size(array);
if s, else, warning('array must be an array'), end    %warning if s is empty
if idx
    q = array(idx);
    qshift = circshift(q,[0,-1]);
    indices = find(q~=qshift);
    if indices, ret = idx(indices);, else, ret = length(q);, end
else
    array=array;
    arrayshift = circshift(array,[0,-1]);
    indices = find(array~=arrayshift);
    if indices, ret = indices;, else, ret = length(array);, end
end