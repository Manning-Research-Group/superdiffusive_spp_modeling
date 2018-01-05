function [tolerance_check sibling_cell] = div_param_check(xyzs_id, xyzs_id_columns, event_info, cell_mat, tolerance, radius_check)

if nargin < 7
    radius_check = 40; %pixel radius to look for other cells
end

divided_cell_pos = [cell_mat(1,1) cell_mat(1,2)];

% Let's check cells within some distance and see if they match
bool1 = xyzs_id(:,xyzs_id_columns - 1) == event_info(1,4); %boolean for division frame
bool2 = sqrt((xyzs_id(:,1)-divided_cell_pos(1)).^2+(xyzs_id(:,2)-divided_cell_pos(2)).^2) <= radius_check;
bool3 = xyzs_id(:,xyzs_id_columns) ~= cell_mat(xyzs_id_columns);
xyzs_id_pot_cell = xyzs_id((bool1&bool2&bool3),:); % cells within radius and in frame

% Now let's check to see if any of those cells are within the tolerance for
% area (column 6)
error = zeros(size(xyzs_id_pot_cell,1),1);
for j = 1:size(xyzs_id_pot_cell,1)
    error(j) = abs(cell_mat(6)-xyzs_id_pot_cell(j,6))/mean([cell_mat(6) xyzs_id_pot_cell(j,6)]);
end
if min(error) < tolerance %If cell closest in area is within tolerance, let's assign that cell the sibling cell
   tolerance_check = 1;
   [value, index] = min(error);
   sibling_cell = xyzs_id_pot_cell(index, xyzs_id_columns); 
else
    tolerance_check = 0;
    sibling_cell = 0;
end
