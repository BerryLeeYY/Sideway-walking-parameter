function filtered_velocity = SW_filter_func(velocity_result)
threshold = 50000;
for i = 1:length(velocity_result)
    if abs(velocity_result(i)) < threshold 
        velocity_result(i) = 0;
    else
         velocity_result(i) =  velocity_result(i);
    end
end

filtered_velocity = velocity_result;
end