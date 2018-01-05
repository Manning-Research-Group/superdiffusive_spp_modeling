function [xshort_slope, xlong_slope, yshort_slope, ylong_slope, alpha2x, alpha2y, alpha2b, tMSD, intersect_x, intersect_y, intercept_upper_x, intercept_upper_y ] = MSD_henderson_theta(xyzs_id, theta, new_dir, filename, plot_toggle, frameindx, cellindx, framerate)
% Function that plots cell tracks and then calls on a subfunction to 
% calculate the decomposed MSD profiles for cell tracks. In that 
% subfunction a line is fit to the short and long timescale behavior and
% the intercept of the long timescale behavior is calculated (this is 
% defined as a mobility parameter in our paper).
%
%  2/8/2012
%  M. Brasch, R. Baker, L. Manning
% 
%  INPUTS:
%  xyzs_id: matrix of tracked cell information for all frames after
%           post-processing
%  theta: angle of orientation (wrinkle direction)
%  new_dir: directory for saving figures
%  filename: prefix name for saving the tracks and MSD figures
%  plot_toggle: 1-plotting on; 0-plotting off
%  frameindx: column number containing frames
%  cellindx: column number containing cell IDs
%  framerate: how many minutes between frames
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

if nargin < 5;
    plot_toggle = 1;
end
if nargin < 6;
    frameindx = 12;
end
if nargin < 7;
    cellindx  = 13;
end
if nargin < 8;
    framerate = 3;
end

nframes = max(xyzs_id(:,frameindx));
ncells = max(xyzs_id(:,cellindx));

% initialization
if plot_toggle ==1
    figure(4);
    cmap = colormap(lines(ncells));
end
ltracks = zeros(ncells,1);
xpos = zeros(ncells,nframes);
ypos = zeros(ncells,nframes);
zpos = zeros(ncells,nframes);

% Plot tracks
for i=1:ncells
    %boolean matrix identifying all cells with id i
    boolcell = (xyzs_id(:,cellindx) == (i));
    ltracks(i) = nnz(boolcell); % length of that track
    if plot_toggle ==1
        hold on;
        plot(xyzs_id(boolcell,1),xyzs_id(boolcell,2), '-','Linewidth',2,'Color',cmap(i,:));
    end
    
    xpos(i,xyzs_id(boolcell,frameindx)) = xyzs_id(boolcell,1);
    ypos(i,xyzs_id(boolcell,frameindx)) = xyzs_id(boolcell,2);
    %zpos(i,xyzs_id(boolcell,9)) = xyzs_id(boolcell,3);     
    
end

% Save figure
if plot_toggle == 1;
%     sdf(4,'TrackingPaper');
    save_name4 = [new_dir, '\', filename, '_full_tracks'];
%     export_fig(save_name4, '-eps');
    h = gcf;
    saveas(h,save_name4, 'fig');
end

%time (in minutes) between frames
time = framerate*(1:nframes);

[ xshort_slope, xlong_slope, yshort_slope, ylong_slope, alpha2x, alpha2y, alpha2b, tMSD, intersect_x, intersect_y, intercept_upper_x, intercept_upper_y] = get_MSD_only_tracks_theta3( xpos', ypos', zpos',time, new_dir, filename, theta, plot_toggle);


end

