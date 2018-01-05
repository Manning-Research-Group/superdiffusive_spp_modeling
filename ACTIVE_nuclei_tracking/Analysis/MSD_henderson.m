function [short_slope, long_slope, intersect, intercept_upper] = MSD_henderson(xyzs_id, new_dir, filename, plot_toggle, frameindx, cellindx,framerate)
% Function that plots cell tracks and then calls on a subfunction to 
% calculate the MSD profiles for cell tracks. In that subfunction a line is
% fit to the short and long timescale behavior and the intercept of the 
% long timescale behavior is calculated (this is defined as a 
% mobility parameter in our paper).
%
%  2/8/2012
%  M. Brasch, R. Baker, L. Manning
% 
%  INPUTS:
%  xyzs_id: matrix of tracked cell information for all frames after
%           post-processing
%  new_dir: directory for saving figures
%  filename: prefix name for saving the tracks and MSD figures
%  plot_toggle: 1-plotting on; 0-plotting off
%  frameindx: column number containing frames
%  cellindx: column number containing cell IDs
%  framerate: how many minutes between frames
%
%  OUTPUTS:
%  short_slope: short timescale slope of logMSD-logt data
%  long_slope: long timescale slope of logMSD-logt data
%  intersect: coordinates for the intersect of short and long timescale
%             slopes
%  intercept_upper: mobility parameter defined as the y-intercept of the
%                   line fit to the long timescale data
%
% 
if nargin < 4
    plot_toggle = 1;
end
if nargin < 5
    frameindx = 12;
end
if nargin < 6
    cellindx = 13;
end
if nargin < 13
    framerate = 3;
end

nframes = max(xyzs_id(:,frameindx));
ncells = max(xyzs_id(:,cellindx));

% initialization
if plot_toggle == 1
    figure(1)
    cmap = colormap(lines(ncells));
end
ltracks = zeros(ncells,1);
xpos = zeros(ncells,nframes);
ypos = zeros(ncells,nframes);
zpos = zeros(ncells,nframes);
%tstartfin = zeros(ncells,2);


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

% Plot Tracks
if plot_toggle == 1
%     sdf(1,'TrackingPaper');
    save_name1 = [new_dir, '\', filename, '_full_tracks'];
%     export_fig(save_name1, '-eps');
    h = gcf;
    saveas(h,save_name1, 'fig');
end

%time (in minutes) between frames
time = framerate*(1:nframes);

[short_slope long_slope MSD, NMSD, alpha2, tMSD, intersect, intercept_upper ] = get_MSD_only_tracks2( xpos', ypos', zpos', new_dir, filename,time,plot_toggle);



end

