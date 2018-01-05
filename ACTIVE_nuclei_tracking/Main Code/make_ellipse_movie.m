function [a ] = make_ellipse_movie(xyzs_id, xyzs_id_columns, stackname, data_folder, figbase, folder_name)
%Function to generate a movie overlay of cell images and ellipse
%information. Sibling information is plotted in different colors as cells
%divide. 
%
%  7/7/2012
%  R. Baker, M. Brasch, L. Manning
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
%  xyzs_id: matrix of particle information post tracking and correcting
%           (output from collision_tags)
%  filename: name of the image stack used
%  figbase: prefix to call the figure
%  folder_name: name of folder where images will be saved
%
%  OUTPUTS:
%  a: let user know function is complete
%
%  LOCAL PARAMETERS:


if nargin < 6
    folder_name = 'Overlay_Images';
end
if nargin < 5
    figbase = 'fig';
end

xc = xyzs_id(:,1);      %x center of mass
yc = xyzs_id(:,2);      %y center of mass
a = xyzs_id(:,3);       %major axis
b = xyzs_id(:,4);       %minor axis
theta= xyzs_id(:,5);    %angle
npoints = 25;           %number of points to generate (to plot ellipse)
pel = cell(length(xyzs_id),1);

%Generate ellipse data from xc, yc, a, b, and theta. npoints varies
%depending on user preference
for i = 1:length(xyzs_id)
    eplotdata = plot_ellipse(xc(i),yc(i),a(i),b(i),theta(i),npoints);
    pel{i} = eplotdata;
end

image_info = imfinfo(stackname);
number_images = numel(image_info);

mkdir(data_folder, folder_name) % make folder to store video images
type_chk = strfind(stackname, '\');
if ~isempty(type_chk)
    addpath([data_folder, '\', folder_name]);
else
    addpath([data_folder, '/', folder_name]);
end

% plot frame by frame
for j = 1:number_images

    A = imadjust(imread(stackname, j, 'Info', image_info));
    
    h = figure(j);
    hold off;
    
    % Plot grayscale image
    imagesc(A)
    colormap('gray')
    hold on
    
    % Find all particles for each frame
    boolpid = xyzs_id(:,xyzs_id_columns - 1) == j;        
    
    % Number of particles in this frame
    num_p_frame = length(find(boolpid==1));     

    bool_list = find(boolpid==1);
    
    %%%% Plot all particles for frame k
    for k=1:num_p_frame
        es=pel{bool_list(k)};
        plot(es(1,:),es(2,:),'-r', 'LineWidth', 2)
        
        % Plot siblings as separate particles (y,m,c,g,b,k alternating); 
        % single particles appear red.
%         switch xyzs_id(bool_list(k),xyzs_id_columns + 2)
%             case 1
%                 plot(es(1,:),es(2,:),'-y') 
%             case 2
%                 plot(es(1,:),es(2,:),'-m') 
%             case 3
%                 plot(es(1,:),es(2,:),'-c') 
%             case 4
%                 plot(es(1,:),es(2,:),'-g')
%             case 5
%                 plot(es(1,:),es(2,:),'-b')
%             case 6
%                 plot(es(1,:),es(2,:),'-k')
%             otherwise
%                 plot(es(1,:),es(2,:),'-r') 
%         end
    end
    
    figname = [figbase, '_',num2str(j)];
    %figname = strcat('fig_',num2str(j));
    if ~isempty(type_chk)
        saveas(h, [data_folder '\' folder_name '\' figname] , 'tif');
    else
        saveas(h, [data_folder '/' folder_name '/' figname] , 'tif');
    end
    close 

end

a='done';

end

 