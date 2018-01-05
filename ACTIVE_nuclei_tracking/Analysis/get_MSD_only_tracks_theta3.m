function [xshort_slope, xlong_slope, yshort_slope, ylong_slope, alpha2x, alpha2y, alpha2b, tMSD, intersect_x, intersect_y, intercept_upper_x, intercept_upper_y ] = get_MSD_only_tracks_theta3( xpos, ypos, zpos, time, new_dir, filename, theta, plot_toggle)
% Function that calculates the MSD profiles for cell tracks.  The MSD 
% profiles are decomposed into x and y directions to account for 
% anisotropy. A line is fit to the short and long timescale behavior and 
% the intercept of the long timescale behavior is calculated (this is 
% defined as a mobility parameter in our paper).
%
%  5/1/2013
%  M. Brasch, R. Baker, L. Manning
% 
%  INPUTS:
%  xpos,ypos: vector data of x and y coordinates
%  new_dir: directory to save images to
%  filename: name prefix for saving files
%  time: vector data of time points
%  plot_toggle: 1-plotting on, 0-plotting off
%
%  OUTPUTS:
%  xshort_slope,yshort_slope: short timescale slope of logMSD-logt data for
%                             x and y directions
%  xlong_slope, ylong_slope: long timescale slope of logMSD-logt data for x
%                             and y directions
%  alpha2x, alphs2y, alpha2b: non-gaussian distribution parameter for x, y
%                             and combined directions
%  tMSD: time vector for MSD plotting
%  intersect_x,intersect_y: coordinates for the intersect of short and long 
%                           timescale slopes for x and y directions
%  intercept_upper_x, intercept_upper_y: mobility parameter defined as the 
%                           y-intercept of the line fit to the long 
%                           timescale data
%
% LOCAL PARAMETERS:
throw_away = 60; % Number of end points to disregard due to high stdev
startlow = 1; % Point to start analysis of short timescale slope
starthigh = 5; % Point to finish analysis of short timescale slope

if nargin < 7
    theta = 0;
end
if nargin < 8
    plot_toggle = 1;
end

[Nt Ncells] = size(xpos);
% Now calculate the mean squared displacement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
Nsub = round(14*Nt/15);
dt = time(1);
        tstart = max(Nt-Nsub,1);
        delt = Nt-tstart;
        MSDx = zeros(delt,1);
        M4Dx = zeros(delt,1);
        MSDy = zeros(delt,1);
        M4Dy = zeros(delt,1);
        M4Db = zeros(delt,1);
        MSDb = zeros(delt,1);
        % Mean squared displacement:
        tcounter = 0;
        tlist = 1:delt+1;
        NMSD = zeros(delt+1,1);
        for t1=tstart:(Nt-1)
            for t2=t1+1:Nt
                mdt = t2-t1; % delta t (time difference between t1&t2)
                
                boolt1 = xpos(t1,:) ~= 0;
                boolt2 = xpos(t2,:) ~= 0;
                boolboth = boolt1 & boolt2;% don't include data where positions are zero
                                            % (this happens a lot in tracking data from experiments)
                            % Now we want to project the x,y, components of the velocity
            % vector so that Mcx is in direction of grooves (theta)
            newx1 = xpos(t1,boolboth)*cos(theta) +  ypos(t1,boolboth)*sin(theta);
            newy1 = -xpos(t1,boolboth)*sin(theta) +  ypos(t1,boolboth)*cos(theta);
            newx2 = xpos(t2,boolboth)*cos(theta) +  ypos(t2,boolboth)*sin(theta);
            newy2 = -xpos(t2,boolboth)*sin(theta) +  ypos(t2,boolboth)*cos(theta);
                                        
                % total net displacement between times t1 and t2
            Dispsqx = (newx1 -newx2).^2;
            Dispsqy = (newy1 -newy2).^2;
            Dispboth = Dispsqx + Dispsqy;
            
                % total net displacement between times t1 and t2
                %Disp =  (xpos(t1,boolboth) - xpos(t2,boolboth)).^2 + (ypos(t1,boolboth) - ypos(t2,boolboth)).^2 + (zpos(t1,boolboth) - zpos(t2,boolboth)).^2;
                dispfx = Dispsqx(isfinite(Dispsqx));
                dispfy = Dispsqy(isfinite(Dispsqy));
                dispfb = Dispboth(isfinite(Dispboth));
                if ~isempty(dispfx)
                    % use vectorization to speed up Matlab.
                    % MSD calculates the average net squared distance between all the particles
                    % i.e. sum(dispf)/length(dispf), and adds it, properly
                    % weighted, to the running average
                    MSDx(mdt) = (NMSD(mdt)*MSDx(mdt) + sum(dispfx))/(NMSD(mdt) +length(dispfx));
                    M4Dx(mdt) = (NMSD(mdt)*M4Dx(mdt) + sum(dispfx.^2))/(NMSD(mdt) +length(dispfx));
                    MSDy(mdt) = (NMSD(mdt)*MSDy(mdt) + sum(dispfy))/(NMSD(mdt) +length(dispfy));
                    M4Dy(mdt) = (NMSD(mdt)*M4Dy(mdt) + sum(dispfy.^2))/(NMSD(mdt) +length(dispfy));
                    
                    MSDb(mdt) = (NMSD(mdt)*MSDb(mdt) + sum(dispfb))/(NMSD(mdt) +length(dispfb));
                    M4Db(mdt) = (NMSD(mdt)*M4Db(mdt) + sum(dispfb.^2))/(NMSD(mdt) +length(dispfb));
                    
                    NMSD(mdt) = NMSD(mdt)+length(dispfx);
                end
            end
        end
    
    for p = 1:length(NMSD)-1;
        stdMSD(p) = 1/sqrt(NMSD(p));
    end
        
    %"non-Gaussian" parameter    
    alpha2x = (3/5)*(M4Dx./(MSDx.^2)) -1;
    alpha2y = (3/5)*(M4Dy./(MSDy.^2)) -1;
    alpha2b = (3/5)*(M4Db./(MSDb.^2)) -1;
          
    % Ravg is a characteristic length scale (for experimental data =1 pixel)
    Ravg = 1;
    Nkeep = length(MSDx);
    tMSD = time(2:Nkeep+1);
    
    
    % plot the "non-Gaussian" parameter to look for cage breaking
    if plot_toggle == 1;
        h = figure(5);
        semilogx(tMSD, alpha2x,'r-o');
        hold on;
        semilogx(tMSD, alpha2y,'b-o');
        semilogx(tMSD, alpha2b, 'm-o');
        hold off

%         sdf(5,'TrackingPaper');
        save_name5 = [new_dir, '\', filename, '_alpha_no_fit'];
%         export_fig(save_name5, '-eps');
        saveas(h,save_name5, 'fig');
    end
    
    coeff1 = [1; -1]; % Starting guess for x(1), x(2)
    lb = [0 -Inf]; 
    ub = [Inf 0];
    [x,resnorm] = lsqcurvefit(@mypowerfun,coeff1,tMSD,alpha2x', lb, ub);
    coeff2 = [1; -1]; % Starting guess for x(1), x(2)
    [x2,resnorm2] = lsqcurvefit(@mypowerfun,coeff2,tMSD,alpha2y', lb, ub);

    if plot_toggle == 1;
        h = figure(6);
        semilogx(tMSD, alpha2x,'r-o');
        hold on
        semilogx(tMSD, alpha2y,'b-o');
        plot(tMSD,mypowerfun(x,tMSD), 'g', 'Linewidth', 2);
        strunits = get(h, 'Unit');
        str1 = ['x parameter a is: ', num2str(x(1))];
        text(0.3,0.8, str1,'Units','normalized');
        set(h, 'Unit', strunits);
        str2 = ['x parameter b is: ', num2str(x(2))];
        text(0.3,0.75, str2,'Units','normalized');

        semilogx(tMSD, alpha2y,'b--o');
        plot(tMSD,mypowerfun(x2,tMSD), 'k', 'Linewidth', 2)
        str1 = ['y parameter a is: ', num2str(x2(1))];
        text(0.3,0.65, str1,'Units','normalized');
        set(h, 'Unit', strunits);
        str2 = ['y parameter b is: ', num2str(x2(2))];
        text(0.3,0.60, str2,'Units','normalized');

        xlabel('time');
        ylabel('\alpha_2');

%         sdf(6,'TrackingPaper');
        save_name6= [new_dir, '\', filename, '_alpha_power_fit'];
%         export_fig(save_name6, '-eps');    
        saveas(h,save_name6, 'fig');
    end

    if max(MSDx ~= 0) %i.e. there's no problem with the data
        
        for v=1:length(MSDx)
            errorlogx(v) = 1/sqrt(NMSD(v))*log10(MSDx(v));
            percent_errorx(v) = 100*errorlogx(v)/log10(MSDx(v));
            errorlogy(v) = 1/sqrt(NMSD(v))*log10(MSDy(v));
            percent_errory(v) = 100*errorlogy(v)/log10(MSDy(v));
        end
        
        

        % this characterizes "upper part" and "lower part"
        % of the MSD in case there is a change in slope with timescale
        finallow = round(0.4*(Nkeep-throw_away));
        finalhigh = Nkeep-throw_away;
        nlow = finallow;
        nhigh = finalhigh;
        [p1 S1] = polyfit(log10((finallow:finalhigh)*dt),log10(MSDx(finallow:finalhigh)'/(Ravg.^2)),1);
        ntime1 = log10((finallow:finalhigh)*dt);
        upm = p1(1);
        upb = p1(2);
        xlong_slope = p1(1);

        nlow2 = startlow;
        nhigh2 = starthigh;
        [p2 S2] = polyfit(log10((startlow:starthigh)*dt),log10(MSDx(startlow:starthigh)'/(Ravg.^2)),1);
        ntime2 = log10((startlow:starthigh)*dt);
        lowm = p2(1);
        lowb = p2(2);
        xshort_slope = p2(1);
        
        Ax = [1 -upm; 1 -lowm]; b=[upb lowb];
        intersect_x = Ax/b;
        intercept_upper_x = upb;

        [p3 S3] = polyfit(log10((finallow:finalhigh)*dt),log10(MSDy(finallow:finalhigh)'/(Ravg.^2)),1);
        ntime3 = log10((finallow:finalhigh)*dt);
        upm = p3(1);
        upb = p3(2);
        ylong_slope = p3(1);

        nlow = startlow;
        nhigh = starthigh;
        [p4 S4] = polyfit(log10((startlow:starthigh)*dt),log10(MSDy(startlow:starthigh)'/(Ravg.^2)),1);
        ntime4 = log10((startlow:starthigh)*dt);
        lowm = p4(1);
        lowb = p4(2);
        yshort_slope = p4(1);
        
        Ay = [1 -upm; 1 -lowm]; b=[upb lowb];
        intersect_y = Ay/b;
        intercept_upper_y = upb;

        % Best MSD plot showing the best fit slopes + diffusion constant
        if plot_toggle == 1;
            h= figure(7);
            hold off;
            errorbar(log10((1:delt)*dt),log10(MSDx'/(Ravg.^2)),errorlogx,'b.')
            hold on;
            errorbar(log10((1:delt)*dt),log10(MSDy'/(Ravg.^2)),errorlogy,'m.')

            %fit for x and y
            plot(ntime1, p1(1)*ntime1 + p1(2), 'k --o') 
            str2x = [ 'Long timescale x slope: ', num2str(p1(1))];

            strunits = get(h, 'Unit');
            set(h, 'Unit', 'normalized');
            text(0.05,0.9, str2x, 'Units','normalized');
            set(h, 'Unit', strunits);
            hold on;
            plot(ntime2, p2(1)*ntime2 + p2(2), 'g --o') 
            strx = [ 'Short timescale x slope: ', num2str(p2(1))];

            strunits = get(h, 'Unit');
            set(h, 'Unit', 'normalized')
            text(0.05,0.8, strx,'Units','normalized');
            set(h, 'Unit', strunits);  

            plot(ntime3, p3(1)*ntime3 + p3(2), 'k --o') 
            str2y = [ 'Long timescale y slope: ', num2str(p3(1))];

            strunits = get(h, 'Unit');
            set(h, 'Unit', 'normalized');
            text(0.05,0.7, str2y, 'Units','normalized');
            set(h, 'Unit', strunits);
            hold on;
            plot(ntime4, p4(1)*ntime4 + p4(2), 'g --o') 
            stry = [ 'Short timescale y slope: ', num2str(p4(1))];

            strunits = get(h, 'Unit');
            set(h, 'Unit', 'normalized')
            text(0.05,0.6, stry,'Units','normalized');
            set(h, 'Unit', strunits);  

            xlabel( 'log_{10}(time) (mins)')
            ylabel( 'log_{10}(distance) (microns)')
            axis([0 3.5 0 6])

%             sdf(7,'TrackingPaper');
            save_name7= [new_dir, '\', filename, '_MSD_decomp'];
%             export_fig(save_name7, '-eps');
            saveas(h,save_name7, 'fig');
        end

    else
        
        badMSD = 'uh oh, MSD is zero for whole time';
        diffconst = -1;
        upm = -1;
        lowm = -1;
        intx = -1;
        inty = -1;
    end
    
end

