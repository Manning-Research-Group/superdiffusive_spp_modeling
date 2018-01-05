function [done] = gyration_plot2(xyzs_id,xyzs_id_columns, gyration_angle,new_dir,filename)
% Function that plots cell postions, rotated by the principle axis of the
% largest eighenvector, with respect to a common origin. A plot of the
% rotated terminal positions is also generated.
%
%  6/25/2013
%  R. Baker, M. Brasch
%
%  INPUTS:
%  xyzs_id: matrix of tracked cell information for all frames after
%           post-processing
%  xyzs_id_columns: column number containing cell IDs
%  gyration_angle: vector of principle axis (largest eigenvector) for each
%           cell
%  new_dir: directory for saving figures
%  filename: prefix name for saving plots
%
%  OUTPUTS: None
%
% Sort matrix by cell ID and frame
xyzs_id = sortrows(xyzs_id,[xyzs_id_columns xyzs_id_columns-1]);
cell_IDs = unique(xyzs_id(:,xyzs_id_columns));

% Preallocation
frame_cell = cell(max(xyzs_id(:,xyzs_id_columns-1))-1,1);
index = ones(max(xyzs_id(:,xyzs_id_columns-1)),1);

% Preallocate final position information for each cell
end_plot = zeros(size(cell_IDs,1),2);

% Make sure all angles are from 0 to pi
gyration_angle(gyration_angle(:) < 0) = gyration_angle(gyration_angle(:) < 0) + pi; 

% Calculate rotated positions for each cell in every frame
for i=1:size(cell_IDs,1);
    % Isolate cell information
    bool = xyzs_id(:,xyzs_id_columns) == cell_IDs(i);
    cell_mat = xyzs_id(bool,:);
    xpos = cell_mat(:,1);
    ypos = cell_mat(:,2);
    % Rotate coordinates by the principle axis
    cell_mat(:,1) = xpos*cos(gyration_angle(i)) + ypos*sin(gyration_angle(i));
    cell_mat(:,2) = -xpos*sin(gyration_angle(i)) + ypos*cos(gyration_angle(i));
    % Give cell a common origin of 0,0
    cell_mat(:,1) = cell_mat(:,1) - cell_mat(1,1);
    cell_mat(:,2) = cell_mat(:,2) - cell_mat(1,2);
    
    if size(cell_mat,1) < 2
        continue
    end
    
    % Determine cell positions for every time step and store
    for k=2:size(cell_mat,1);
        delt_f = cell_mat(k,xyzs_id_columns-1) - cell_mat(1,xyzs_id_columns-1);
        frame_cell{delt_f}(index(delt_f),1) = cell_mat(k,1);
        frame_cell{delt_f}(index(delt_f),2) = cell_mat(k,2);
        index(delt_f) = index(delt_f)+1;
    end
    
    % Store terminal position information
    end_plot(i,1) = cell_mat(size(cell_mat,1),1);
    end_plot(i,2) = cell_mat(size(cell_mat,1),2);

end

new_dir2 = [new_dir, '\', filename, '_gyration_plots'];
mkdir(new_dir2);
% Generate figure of all cell positions for each timestep (h)
for h=1:size(frame_cell,1);
    if isempty(frame_cell{h})
        continue
    end
    
    b = figure(9);
    scatter(frame_cell{h}(:,1),frame_cell{h}(:,2), 10);
    axis([-1000 1000 -1000 1000])
    title(['Gyration Plot Frame ',num2str(h)])
    xlabel('Largest Eigenvector Direction')
    ylabel('Smallest Eigenvector Direction')    
%     sdf(9,'TrackingPaper');
    save_name9 = [new_dir2, '\', filename, '_gyration_plot', num2str(h)];
%     export_fig(save_name9, '-tif');
    saveas(b,save_name9, 'tif');
end

% Generate figure of terminal cell positions
n = figure(10);
scatter(end_plot(:,1), end_plot(:,2));
axis ([-1000 1000 -1000 1000]);
title('Final Cell Location (Post Principle Axis Rotation)')
xlabel('Largest Eigenvector Direction')
ylabel('Smallest Eigenvector Direction')

% sdf(10,'TrackingPaper');
save_name10 = [new_dir, '\', filename, '_gyration_plot'];
% export_fig(save_name10, '-eps');
saveas(n,save_name10, 'fig');

done=1;