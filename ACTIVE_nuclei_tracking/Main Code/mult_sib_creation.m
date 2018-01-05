function [mult_sib_array, chk_mult_sib] = mult_sib_creation(xyzs_id,xyzs_id_columns)

index = 1;

for i = 1:max(xyzs_id(:,xyzs_id_columns-1));
    bool = xyzs_id(:,xyzs_id_columns-1) == i;
    if sum(bool) > 0;
        event_num = max(xyzs_id(bool,xyzs_id_columns - 5)); 
        
        for j=1:event_num
            bool2 = xyzs_id(:,xyzs_id_columns-1) == i & xyzs_id(:,xyzs_id_columns-5) == j;
            mult_sib_array{index,1} = i;
            mult_sib_array{index,2} = j;
            mult_sib_array{index,3} = xyzs_id(bool2,xyzs_id_columns);
            mult_sib_array{index,4} = 0;
            index = index + 1;
        end
    end
end

%check the mult_sib_array exists:
chk_mult_sib = exist('mult_sib_array','var');
if chk_mult_sib == 0
   mult_sib_array = []; 
end