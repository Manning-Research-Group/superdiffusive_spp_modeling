function [cell_count_first, cell_count_last, total_cells] = num_cells(xyzs_id, xyzs_id_columns, n_divisions)
    
    %Total number of cells in first frame
    bool_first = xyzs_id(:,xyzs_id_columns-1) == 1;
    cell_count_first = sum(bool_first)-n_divisions;

    %Total number of cells in final frame
    bool_last = xyzs_id(:,xyzs_id_columns-1) == 480;
    cell_count_last = sum(bool_last);
    
    %Find the total number of cells across all frames
    total_cells = size(unique(xyzs_id(:,xyzs_id_columns)),1);
    
end