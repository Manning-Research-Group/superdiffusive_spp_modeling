function [ xa xtime xc ya ytime yc Mcx Mcy NMSD timewindows fstart xv yv] = velocity_auto_theta( xyzs_id, xyzs_id_columns, new_dir, filename, theta, framerate, plot_toggle, vx_col, vy_col)
%Function to calculation velocity autocorrelation for tracked cell data
% M. Lisa Manning 4 March 2013
%  INPUTS: xyzs_id - matrix with tracked cell information from
%                    run_tracking_contour
%          xyzs_id_columns - scalar variable that is the column for cell
%                           ids (frame column is assumed to be one column to
%                           the left.
%          framerate - amount of time between frames
%
%
%  OUTPUTS:
%      Mcx: correlation in x
%      Mcy: correlation in y
%      NMSD: number of data points contributing to each entry in Mcx,Mcy
%      timewindows: the delta t associated with each entry in Mcx,Mcy
%      fstart: the starting frame for the calculation
%      xv(yv): matrix of x(y)-components of velocity vector
%

if nargin < 5
    theta = 0;
end

if nargin < 6
    framerate = 3;
end

if nargin < 7
    plot_toggle = 1;
end

if nargin < 8
    vx_col = size(xyzs_id,2)-1;
end

if nargin < 9
    vy_col = size(xyzs_id,2);
end

nframes = max(xyzs_id(:,xyzs_id_columns-1));
ncells = max(xyzs_id(:,xyzs_id_columns));

%Populate xv, yv to be a nframes by ncells matrix of the velocities
xv = zeros(nframes, ncells);
yv = zeros(nframes, ncells);

for i=1:ncells
    boolcell = (xyzs_id(:,xyzs_id_columns) == (i));
    xv(xyzs_id(boolcell,xyzs_id_columns-1),i) = xyzs_id(boolcell,vx_col);
    yv(xyzs_id(boolcell,xyzs_id_columns-1),i) = xyzs_id(boolcell,vy_col);
end

% initialization

% Nsub is the place to start calculating the velocity autocorrelation;
% Skip first data points to get rid of initialization problems
Nsub = round(3*nframes/4);
fstart = max(nframes-Nsub,1);
delf = nframes-fstart;

NMSD = zeros(delf+1,1); %number of entries in each bin
Mcx = zeros(delf+1,1); % correlation in x
Mcy = zeros(delf+1,1); % correlation in y
%Mcxn = zeros(delf+1,1);
%Mcyn = zeros(delf+1,1);
    for f1=fstart:nframes
        for f2=f1+1:nframes
            mdt = f2-f1; % delta f (frame difference between f1&f2)

            boolt1 = xv(f1,:) ~= 0;
            boolt2 = xv(f2,:) ~= 0;
            boolboth = boolt1 & boolt2;% don't include data where velocities are zero
                                        % (this happens a lot in tracking data from experiments)
            % Now we want to project the x,y, components of the velocity
            % vector so that Mcx is in direction of grooves (theta)
            newxv1 = xv(f1,boolboth)*cos(theta) +  yv(f1,boolboth)*sin(theta);
            newyv1 = -xv(f1,boolboth)*sin(theta) +  yv(f1,boolboth)*cos(theta);
            newxv2 = xv(f2,boolboth)*cos(theta) +  yv(f2,boolboth)*sin(theta);
            newyv2 = -xv(f2,boolboth)*sin(theta) +  yv(f2,boolboth)*cos(theta);
                                        
            % total net displacement between times t1 and t2
            cx = newxv1.*newxv2;
            cy = newyv1.*newyv2;
            %cxnorm = cx./sqrt(abs(xv(f1,boolboth)).*abs(xv(f2,boolboth))); 
            %cynorm = cy./sqrt(abs(yv(f1,boolboth)).*abs(yv(f2,boolboth)));
            cxf = cx(isfinite(cx));
            cyf = cy(isfinite(cy));
            %cxfn = cx(isfinite(cxnorm));
            %cyfn = cy(isfinite(cynorm));
            if ~isempty(cxf)
                % use vectorization to speed up Matlab.
                % Mcx calculates the average net velocity correlation between all the particles
                % i.e. sum(dispf)/length(dispf), and adds it, properly
                % weighted, to the running average
                Mcx(mdt) = (NMSD(mdt)*Mcx(mdt) + sum(cxf))/(NMSD(mdt) +length(cxf));
                Mcy(mdt) = (NMSD(mdt)*Mcy(mdt) + sum(cyf))/(NMSD(mdt) +length(cyf));
                %Mcxn(mdt) = (NMSD(mdt)*Mcxn(mdt) + sum(cxfn))/(NMSD(mdt) +length(cxfn));
                %Mcyn(mdt) = (NMSD(mdt)*Mcyn(mdt) + sum(cyfn))/(NMSD(mdt) +length(cyfn));
                NMSD(mdt) = NMSD(mdt)+length(cxf);
            end
        end
    end

% now for some plots    
time = framerate*(1:nframes);
timewindows = framerate*(1:delf+1);
dt = time(1);

lb = [0 -Inf 0]; 
ub = [Inf 0 Inf];
timewindows_fit = timewindows(6:end); %Remove the first few points from the fit
Mcx_fit = Mcx(6:end);   %Remove the first few points from the fit
start_val = max(Mcx_fit);
coeff1 = [start_val; -0.01; 0]; % Starting guess for x(1), x(2), and x(3)
[x,resnorm] = lsqcurvefit(@myexpfun,coeff1,timewindows_fit,Mcx_fit',lb, ub);
xa = x(1); xtime = -x(2); xc = x(3);

Mcy_fit = Mcy(6:end);   %Remove the first few points from the fit
start_val2 = max(Mcy_fit);
coeff2 = [start_val2; -0.01; 0]; % Starting guess for x(1) and x(2)
[x2,resnorm2] = lsqcurvefit(@myexpfun,coeff2,timewindows_fit,Mcy_fit',lb, ub);
ya = x2(1); ytime = -x2(2); yc = x2(3);

if plot_toggle == 1;
    figure(8);
    h = gcf;
    plot(timewindows, Mcx, 'b--o');
    % title('x-velocity data')
    hold on;
    plot(timewindows_fit,myexpfun(x,timewindows_fit), 'g', 'Linewidth', 2)
    strunits = get(h, 'Unit');
    str1 = ['x-vel x(1) is: ', num2str(x(1))];
    text(0.15,0.8, str1,'Units','normalized');
    set(h, 'Unit', strunits);
    str2 = ['x-vel x(2) is: ', num2str(x(2))];
    text(0.15,0.75, str2,'Units','normalized');
    str3 = ['x-vel x(3) is: ', num2str(x(3))];
    text(0.15,0.7, str3,'Units','normalized');
    % hold off;
    % figure(71)
    % title('y-velocity data')
    plot(timewindows, Mcy, 'm--x');
    % hold on;
    plot(timewindows_fit,myexpfun(x2,timewindows_fit), 'k', 'Linewidth', 2)
    % hold off;
    str1 = ['y-vel x(1) is: ', num2str(x2(1))];
    text(0.15,0.65, str1,'Units','normalized');
    set(h, 'Unit', strunits);
    str2 = ['y-vel x(2) is: ', num2str(x2(2))];
    text(0.15,0.60, str2,'Units','normalized');
    str3 = ['y-vel x(3) is: ', num2str(x2(3))];
    text(0.15,0.55, str3,'Units','normalized');
    xlabel('\Delta t');
    ylabel('<v(0)v(\Delta t)>');
    legend('x direction','x exp fit', 'y direction', 'y exp fit');
    axis([0 1200 -0.05 1.05])
    hold off;

%     sdf(8,'TrackingPaper');
    save_name8= [new_dir, '\', filename, '_velocity'];
%     export_fig(save_name8, '-eps');
    saveas(h,save_name8, 'fig');
end

end

