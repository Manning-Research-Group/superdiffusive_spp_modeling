function s = bpass(img,lambda,w)
%      7-23-03  Maria Kilfoil
%
% 		Implements a real-space bandpass filter to suppress pixel noise and 
%       slow-scale image variations while retaining information of a characteristic size.
% 		*Works with anisotropic 3d cube data*
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
%  CALLING SEQUENCE:
% 		res = bpass(img, lambda, w)
%  INPUTS:
% 		img:	two-dimensional array to be filtered.
% 		lambda: characteristic lengthscale of noise in pixels. Additive noise averaged 
%                   over this length should vanish. May assume any positive floating value.
% 			Make it a 3-vector if aspect ratio is not 1:1:1.
% 		w: A length in pixels somewhat larger than *half* a typical object. Must be an odd valued 
%                   integer. Make it a 3-vector if aspect ratio is not 1:1:1.
%  OUTPUTS:
% 		res:	filtered image.
%  PROCEDURE:
% 		simple 'mexican hat' wavelet convolution yields spatial bandpass filtering.
%  NOTES:
%		based on "Methods of digital video microscopy for colloidal studies", John Crocker 
%       and David Grier, J. Colloid Interface Sci. 179, 298 (1996), and on bpass.pro IDL code 
%       written by John Crocker and David Grier. 
%
clear s t
a=double(img);
b=double(lambda);
w = round(max(w,2*b));
N = 2*w + 1;
r = [-w:w]/(2*b);
xpt = exp(-r.^2);  
B=(sum(xpt))^2;
xpt = xpt / sum(xpt);    %xpt=exp(-[-w:w].^2/4*lambda^2)/sqrt(B)
factor=((sum(xpt.^2))^2-1/N^2)*B;  %sum(xpt.^2)=1/Bexp(-(i^2+j^2)/(4*lambda^2))
% note: N not N^2 etc since doing 2D conv along each axis separately
gx = xpt;    
gy=gx';
bx = zeros(1,N)-1./N;  
by=bx';
% g = conv2( a, gx );
% g = conv2( g, gy );
% b = conv2( a, bx );
% b = conv2( b, by );
g = conv2(gx,gy,a,'valid');
b = conv2(bx,by,a,'valid');
res = g-b;
s=max(res/factor,0);
% whos s
tmp=zeros(length(a(:,1)),length(a(1,:)));
tmp(w+1:(length(s(:,1))+w),w+1:(length(s(1,:))+w))=s;
s=tmp; %filtered image with w width zero border
