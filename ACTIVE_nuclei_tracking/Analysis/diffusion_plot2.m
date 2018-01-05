function [done] = diffusion_plot2(xyzs_id,xyzs_id_columns, theta,new_dir,filename)
% Function to generate diffusion plots of cell locations from a common 
% origin as a function of time. Also generates a plot showing final cell 
% locations from a common origin.
%
%  7/1/2013
%  R. Baker, M. Brasch
%
%  INPUTS:
%  xyzs_id: matrix of tracked cell information for all frames after
%           post-processing
%  xyzs_id_columns: column number containing cell IDs
%  theta: angle of wrinkle direction
%  new_dir: directory to save file to
%  filename: name of the file stack used for saving images
%
%  OUTPUTS: None

% Sort the matrix first by cell number then frame
xyzs_id = sortrows(xyzs_id,[xyzs_id_columns xyzs_id_columns-1]);

% Create vector of unique cell IDs
cell_IDs = unique(xyzs_id(:,xyzs_id_columns));

% Rotate cell positions by theta
xpos = xyzs_id(:,1)*cos(theta) +  xyzs_id(:,2)*sin(theta);
ypos = -xyzs_id(:,1)*sin(theta) +  xyzs_id(:,2)*cos(theta);

xyzs_id(:,1) = xpos;
xyzs_id(:,2) = ypos;

frame_cell = cell(max(xyzs_id(:,xyzs_id_columns-1))-1,1);
index = ones(max(xyzs_id(:,xyzs_id_columns-1)),1);

% Loop that calculates cell positions for every frame a cell is present
for i=1:size(cell_IDs,1);
    bool = xyzs_id(:,xyzs_id_columns) == cell_IDs(i);
    cell_mat = xyzs_id(bool,:);
    cell_mat(:,1) = cell_mat(:,1) - cell_mat(1,1);
    cell_mat(:,2) = cell_mat(:,2) - cell_mat(1,2);
    
    if size(cell_mat,1) < 2
        continue
    end
    
    % Calculate all cell displacements and store in cell array for every
    % delta_f (change in frame)
    for k=2:size(cell_mat,1);
        delt_f = cell_mat(k,xyzs_id_columns-1) - cell_mat(1,xyzs_id_columns-1);
        frame_cell{delt_f}(index(delt_f),1) = cell_mat(k,1);
        frame_cell{delt_f}(index(delt_f),2) = cell_mat(k,2);
        index(delt_f) = index(delt_f)+1;
    end
    
    % Store final cell locations for plot
    diff_plot(i,1) = cell_mat(size(cell_mat,1),1);
    diff_plot(i,2) = cell_mat(size(cell_mat,1),2);

end

% new_dir2 = [new_dir, '\', filename, '_diffusion_plots'];
% mkdir(new_dir2);
% 
% % Generates plot for each frame
% for h=1:size(frame_cell,1);
%     if isempty(frame_cell{h})
%         continue
%     end
%     
%     b=figure(9);
%     scatter(frame_cell{h}(:,1),frame_cell{h}(:,2), 10);
%     axis([-1000 1000 -1000 1000])
%     title(['Diffusion Plot Frame ',num2str(h)])
%     xlabel('X Distance')
%     ylabel('Y Distance')
%     
% %     sdf(9,'TrackingPaper');
%     save_name9 = [new_dir2, '\', filename, '_diffusion_plot', num2str(h)];
%     saveas(b,save_name9,'tif')
% %     export_fig(save_name9, '-tif');
% end

% Generate final cell location plot
figure(10)
scatter(diff_plot(:,1), diff_plot(:,2));
axis ([-1000 1000 -1000 1000]);
title('Diffusion Plot Frame')
xlabel('X Distance')
ylabel('Y Distance')

% sdf(10,'TrackingPaper');
save_name10 = [new_dir, '\', filename, '_diffusion_plot'];
saveas(gcf,save_name10,'fig')
% export_fig(save_name10, '-eps');


done=1;