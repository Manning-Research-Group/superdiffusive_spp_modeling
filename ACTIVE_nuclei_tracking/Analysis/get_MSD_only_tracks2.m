function [short_slope, long_slope, MSD, NMSD, alpha2, tMSD, intersect, intercept_upper ] = get_MSD_only_tracks2( xpos, ypos, zpos, new_dir, filename,time,plot_toggle)
% Function that calculates the MSD profiles for cell tracks. A line is fit
% to the short and long timescale behavior and the intercept of the long
% timescale behavior is calculated (this is defined as a mobility parameter
% in our paper).
%
%  3/18/2013
%  M. Brasch, R. Baker, L. Manning
% 
%  INPUTS:
%  xpos,ypos,zpos: vector data of x,y, and z coordinates
%  new_dir: directory to save images to
%  filename: name prefix for saving files
%  time: vector data of time points
%  plot_toggle: 1-plotting on, 0-plotting off
%
%  OUTPUTS:
%  short_slope: short timescale slope of logMSD-logt data
%  long_slope: long timescale slope of logMSD-logt data
%  MSD: vector of MSD values for each delta time
%  NMSD: number of calculations per delta time
%  alpha2: non-gaussian distribution parameter
%  tMSD: time vector for MSD plotting
%  intersect: coordinates for the intersect of short and long timescale
%             slopes
%  intercept_upper: mobility parameter defined as the y-intercept of the
%                   line fit to the long timescale data
%
% LOCAL PARAMETERS:
throw_away = 60; % Number of end points to disregard due to high stdev
startlow = 1; % Point to start analysis of short timescale slope
starthigh = 5; % Point to finish analysis of short timescale slope

if nargin < 7
    plot_toggle = 1;
end

[Nt Ncells] = size(xpos);
% Now calculate the mean squared displacement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
Nsub = round(14*Nt/15);
dt = time(1);
        tstart = max(Nt-Nsub,1);
        delt = Nt-tstart;
        MSD = zeros(delt,1);
        M4D = zeros(delt,1);
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
                
                % total net displacement between times t1 and t2
                Disp =  (xpos(t1,boolboth) - xpos(t2,boolboth)).^2 + (ypos(t1,boolboth) - ypos(t2,boolboth)).^2 + (zpos(t1,boolboth) - zpos(t2,boolboth)).^2;
                dispf = Disp(isfinite(Disp));
                if ~isempty(dispf)
                    % use vectorization to speed up Matlab.
                    % MSD calculates the average net squared distance between all the particles
                    % i.e. sum(dispf)/length(dispf), and adds it, properly
                    % weighted, to the running average
                    MSD(mdt) = (NMSD(mdt)*MSD(mdt) + sum(dispf))/(NMSD(mdt) +length(dispf));
                    M4D(mdt) = (NMSD(mdt)*M4D(mdt) + sum(dispf.^2))/(NMSD(mdt) +length(dispf));
                    NMSD(mdt) = NMSD(mdt)+length(dispf);
                end
            %ASD(tcounter) = sum(Disp(:,tcounter))/length(Xpini);
            end
        end
    
    for p = 1:length(NMSD)-1;
        stdMSD(p) = 1/sqrt(NMSD(p));
    end
        
    %"non-Gaussian" parameter    
    alpha2 = (3/5)*(M4D./(MSD.^2)) -1;
          
    % Ravg is a characteristic length scale (for experimental data =1 pixel)
    Ravg = 1;
    Nkeep = length(MSD);
    tMSD = time(2:Nkeep+1);
    
    
    % plot the "non-Gaussian" parameter to look for cage breaking
    if plot_toggle == 1;
        h = figure(2);
        semilogx(tMSD, alpha2,'r-o');
        xlabel('time');
        ylabel('\alpha_2');

%         sdf(2,'TrackingPaper');
        save_name2 = [new_dir, '\', filename, '_full_alpha'];
%         export_fig(save_name2, '-eps');
        saveas(h,save_name2, 'fig');
    end
    
    if max(MSD ~= 0) %i.e. there's no problem with the data
        
    for v=1:length(MSD)
        errorlog(v) = 1/sqrt(NMSD(v))*log10(MSD(v));
        percent_error(v) = 100*errorlog(v)/log10(MSD(v));
    end

    % this characterizes "upper part" and "lower part"
    % of the MSD in case there is a change in slope with timescale
    finallow = round(0.4*(Nkeep-throw_away));
    finalhigh = Nkeep-throw_away;
    
    % Best MSD plot showing the best fit slopes
    % Short Timescale
    nlow = finallow;
    nhigh = finalhigh;
    [p1 S1] = polyfit(log10((nlow:nhigh)*dt),log10(MSD(nlow:nhigh)'/(Ravg.^2)),1);
    ntime1 = log10((nlow:nhigh)*dt);
    long_slope = p1(1);
    upm = p1(1);
    upb = p1(2);
    % Long Timescale
    nlow2 = startlow;
    nhigh2 = starthigh;
    [p2 S2] = polyfit(log10((nlow2:nhigh2)*dt),log10(MSD(nlow2:nhigh2)'/(Ravg.^2)),1);
    ntime2 = log10((nlow2:nhigh2)*dt);
    short_slope = p2(1);
    lowm = p2(1);
    lowb = p2(2);
    % Calculate intersect and long_timescale intercept
    A = [1 -upm; 1 -lowm]; b = [upb lowb];
    intersect = A/b;
    intercept_upper = upb;
    
    % Plot MSD profiles with slope fittings
    if plot_toggle == 1;
        h= figure(3);
        hold off;
        errorbar(log10((1:delt)*dt),log10(MSD'/(Ravg.^2)),errorlog,'m.')
        hold on;
        plot(ntime1, p1(1)*ntime1 + p1(2), 'k --o') 
        str2 = [ 'Long timescale slope: ', num2str(p1(1))];

        strunits = get(h, 'Unit');
        set(h, 'Unit', 'normalized');
        text(0.05,0.8, str2, 'Units','normalized');
        set(h, 'Unit', strunits);
        hold on;

        plot(ntime2, p2(1)*ntime2 + p2(2), 'b --o') 
        str = [ 'Short timescale slope: ', num2str(p2(1))];
        strunits = get(h, 'Unit');
        set(h, 'Unit', 'normalized')
        text(0.05,0.7, str,'Units','normalized');
        set(h, 'Unit', strunits);  
        xlabel( 'log_{10}(time) (mins)')
        ylabel( 'log_{10}(distance) (microns)')
        axis([0 3.5 0 6])

%         sdf(3,'TrackingPaper');
        save_name3 = [new_dir, '\', filename, '_full_MSD'];
%         export_fig(save_name3, '-eps');
        saveas(h,save_name3, 'fig');
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

