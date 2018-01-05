function [xc, yc, a, b, theta,es] = fit_ellipse( xvector, yvector, plottoggle )
% Find a least-squares fit of an ellipse to the 2D dataset given in data. 
% Uses the Fitzgibbon, Pilu & Fisher algorithm.
%
%  2/20/2012
%  L. Manning
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
%  xvector: vector of x-ccordinates from contour
%  yvector: vector of y-coordinates from contour
%  plottoggle: 1-plotting on, 0-off
%
%  OUTPUTS:
%  xc,yc: centroid x,y coordinates
%  a,b: major, minor axes
%  theta: angle of orientation
%  es: fit data from pars
%

 if nargin < 3
     plottoggle = 0;
 end
 
 % center the data on the origin
 s = length(xvector);
 
 mx = mean(xvector);
 my = mean(yvector);
 
 cvert = zeros(2,s);
 cvert(1,:) = xvector - mx;
 cvert(2,:) = yvector - my;
 
  % Build the design matrix (holds the data).
  % Note it has one data point per row.
  % Columns express parameters a-f from the form
  % a x^2 + b x y + c y^2 + d x + e y + f = 0.
 
  D = zeros(6, s);
  D(1,:) = cvert(1,:).^2;
  D(2,:) = cvert(1,:).*cvert(2,:);
  D(3,:) = cvert(2,:).^2;
  D(4,:) = cvert(1,:);
  D(5,:) = cvert(2,:);
  D(6,:) = ones(s,1);
  D = D';
  
  % Scatter matrix: Transpose(D) * D (use 'usual matrix operator' ##). Gives a 6x6 symmetric matrix
  S = D'*D;
 
  % Constraint matrix, giving 4ac-4b^2=1 as constraint.
  C = zeros(6,6);
  C(1,3) = -2;
  C(2,2) = 1;
  C(3,1) = -2;
  C = C';
  
  % Solve the generalized eigensystem S a = lambda C a.
  [V, D] = eig(S,C);
  
  evals = real(diag(D));  
  evecs = V(:,isfinite(evals));
  evals = evals(isfinite(evals)); 
  
  % Get the negative eigenvalue corresponding to the actual solution.
  % If not found, check if there is a small positive one. If not, report
  % error.
  [val ndx] = min(evals);

  tol = 1;
  if val > tol
      error('This routine did not find a negative eigenvalue');
  end
  
  %Get fit parameters: the eigenvector associated with the found eigenvalue.
  pars = evecs(:,ndx);
    
  % Get values we will need for center and axes.
  denominator1 = (pars(2)/2)^2-pars(1)*pars(3);
  
  numerator = 2*(pars(1)*(pars(5)/2)^2 + pars(3)*(pars(4)/2)^2 + pars(6)*(pars(2)/2)^2 - 2*(pars(2)/2)*(pars(4)/2)*(pars(5)/2) - pars(1)*pars(3)*pars(6));
  denominator2 = sqrt((pars(1)-pars(3))^2 + pars(2)^2);
  
  % Get center
  xc = (pars(3)*(pars(4)/2)-(pars(2)/2)*(pars(5)/2))/denominator1 + mx;
  yc = (pars(1)*(pars(5)/2)-(pars(2)/2)*(pars(4)/2))/denominator1 + my;
      
  % Semimajor and semiminor axis.      
  a = sqrt(numerator/(denominator1*(denominator2-(pars(1)+pars(3)))));
  b = sqrt(numerator/(denominator1*(-denominator2-(pars(1)+pars(3)))));

  % Angle between major axis and horizontal.
  if (pars(2) == 0) 
    if (pars(1) <= pars(3))
        theta = 0; 
    else
        theta = pi/2;
    end
  else 
    if (pars(1) == pars(3))
       theta = 0; 
    else
      if (pars(1) <= pars(3))
          theta = .5*atan(pars(2)/(pars(1)-pars(3)));
      else
          theta = pi/2 + .5*atan(pars(2)/(pars(1)-pars(3)));
      end
    end
  end
   
  if (a < b)  
    temp = b;
    b = a;
    a = temp;
    theta = theta-pi/2;
  end
  
  % get contour for ellipse
  es = plot_ellipse(xc, yc,a,b,theta,50);
  
  if plottoggle == 1
    plot(xvector, yvector, 'bo');
    hold on;
    plot(es(1,:), es(2,:), 'r-');
    axis square
  end
  
end

