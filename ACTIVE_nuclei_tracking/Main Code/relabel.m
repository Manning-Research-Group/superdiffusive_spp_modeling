function [ xyzs_id ] = relabel(compiled_mastertable, xyzs_id, xyzs_id_columns)
%   Relabels collision IDs based on compiled_mastertable
%   Collisions may result in mislabeled cell IDs. After the cost function
%   or position cost function analysis, a compiled_mastertable is formed for all 
%   noted sibling cells. The relabel function isolates this matrix and 
%   reversibly (by frame) updates xyzs_id correspondingly.
% 
%  2/8/2012
%  M. Brasch, R. Baker, L. Manning
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
%  compiled_mastertable: matrix of merging event information (column 9 is whether or
%               not the event needs IDs switched)
%  xyzs_id: matrix of particle information post-tracking (output from
%           trackmem_new.m)
%  xyzs_id_columns: column of xyzs_id matrix that contains the cell ID info
%
%  OUTPUTS:
%  xyzs_id: updated particle information after relabeling the merging
%           events that are tagged as needing switched
%

bool_zero = compiled_mastertable(:,4)~= 0;
new_compiled_mastertable = compiled_mastertable(bool_zero, :);

frames = sort(unique(new_compiled_mastertable(:,2)),'descend');

for j = 1:size(frames,1)
    
    frame_num = frames(j);
    mini_mat = new_compiled_mastertable(new_compiled_mastertable(:,2)==frame_num, :);
    
    ncells = size(mini_mat,1);
    cell_array = cell(ncells,2);
    
    % Index the cell IDs to be relabeled (this is necessary because
    % relabelling cannot be done in series, as duplicates occur)
    for k = 1:ncells
        cell_array{k,1} = xyzs_id(:,xyzs_id_columns) == mini_mat(k,3) & xyzs_id(:,xyzs_id_columns-1)>= mini_mat(k,2);
        cell_array{k,2} = mini_mat(k,4);
    end
    
    % Relabel cell IDs
    for s = 1:ncells
        xyzs_id(cell_array{s,1},xyzs_id_columns) = cell_array{s,2};
    end
    
end

% for i = 1:size(compiled_mastertable,1)
%         %Cell 1 Index Values
%         cell1_indices = xyzs_id(:,xyzs_id_columns) == new_compiled_mastertable(i,3) & xyzs_id(:,xyzs_id_columns-1)>= compiled_mastertable(i,2);
% %         %Cell 2 Index Values
% %         cell2_indices = xyzs_id(:,xyzs_id_columns) == compiled_mastertable(i,4) & xyzs_id(:,xyzs_id_columns-1)>= compiled_mastertable(i,2);
% %         %Relabel Cell 1 & 2 xyzs_id Values
%         xyzs_id(cell1_indices,xyzs_id_columns) = compiled_mastertable(i,4);
% %         xyzs_id(cell2_indices,xyzs_id_columns) = compiled_mastertable(i,3);
%  
% end


end

