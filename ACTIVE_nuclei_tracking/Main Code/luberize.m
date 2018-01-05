% Used for kilfoil tracking

function [newtracks] = luberize(tracks)
%  Maria Kilfoil
% reassigns the unique ID# to 0,1,2,3...
% /presort will sort on ID# first, then reassign
% start will begin with that ID#

% function returns a new track array


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

ndat=length(tracks(1,:));
% if (keyword_set(presort)) then begin
%     newtracks=tracks(*,sort(tracks(ndat,*)))
% endif else begin
%     newtracks=tracks
% endelse

newtracks=tracks;

u=unq((newtracks(:,ndat))',[]);
ntracks=length(u);
u=[0,u];
for i=1:ntracks,  newtracks(u(i)+1:u(i+1),ndat) = i; end

% if (keyword_set(start)) then newtracks(ndat,*)=newtracks(ndat,*)+start

