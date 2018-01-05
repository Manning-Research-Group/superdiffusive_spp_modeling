function [compiled_mastertable] = master_compile(mastertable, merging_output)

% Delete any entries in mastertable that don't result in a switch; then
% reduce master table to only include 1) 1st frame, 2) last frame, 3) ID1,
% 4) ID2
master_bool = mastertable(:,9) == 1;
mastertable_new = mastertable(master_bool,1:4);

% Duplicate information and switch columns 3 and 4 (ID1 and ID2)
mastertable_new_dup(:,[1 2 3 4]) = mastertable_new(:,[1 2 4 3]);
if ~isempty(merging_output)
    merging_bool = merging_output(:,3) ~= merging_output(:,4); 
    merging_output = merging_output(merging_bool,:);
end

% Concatenate both mastertables and the merging_output
compiled_mastertable = [mastertable_new; mastertable_new_dup; merging_output];

compiled_mastertable = sortrows(compiled_mastertable, -2);
