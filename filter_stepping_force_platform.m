function filter_stepping_force_platform_output = filter_stepping_force_platform(stepping_force_platform)
filter_stepping_force_platform_output = zeros(1,length(stepping_force_platform));
for i = 1:length(stepping_force_platform)
    if abs(stepping_force_platform(i)) < 2
        value = 0;
        filter_stepping_force_platform_output(i) = value;
    else
        filter_stepping_force_platform_output(i) = stepping_force_platform(i);
    end
end
end