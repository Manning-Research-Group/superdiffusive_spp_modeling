angle_vec = [0 74.1 74.1 74.1 74.1 0 86.7 86.7 86.7 86.7 0 85.6 85.6 85.6 85.6 0 0 0 0 0 170.7 170.7 170.7 170.7 170.7 88.6 88.6 88.6 88.6 88.6];

angle_vec2 = zeros(1,size(angle_vec,2));
for i = 1:size(angle_vec,2);
    value = angle_vec(i);
    if value == 0;
        angle_vec2(i) = 0;
    else
        angle_vec2(i) = abs(180-value);
    end
end

final_angles = angle_vec2*(pi/180);

