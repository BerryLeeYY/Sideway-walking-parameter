function velocity_result = velocity_func(marker_trajectory)

velocity_result = zeros(length(marker_trajectory));

for i = 1:(length(velocity_result)-1)
    velocity_result(i) = (marker_trajectory(i+1) - marker_trajectory(i))/(1/length(marker_trajectory));
end
end



