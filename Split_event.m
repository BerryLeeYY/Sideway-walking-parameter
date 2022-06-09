%clear;
subject = ('S6_Gait_sw_0001');
subjfile = [(subject),'.mat'];
load(subjfile);
label = Gait_sw_0001.Trajectories.Labeled.Labels;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Right Leg %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Trajectory extraction
RDM5_trajectory = Gait_sw_0001.Trajectories.Labeled.Data(find(strcmp(label, 'RDM5')),1, :);
RDM5_trajectory = reshape(RDM5_trajectory, [length(RDM5_trajectory), 1]);
%%%%%% velocity calculation
marker_trajectory = RDM5_trajectory;
velocity_result = velocity_func(marker_trajectory);
%%%%%% velocity normalized
velocity_result = velocity_result;
normalized_velocity = normalize(velocity_result);
%%%%%% velocity filtering
normalized_velocity = normalized_velocity;
filtered_velocity = filter_func(normalized_velocity);
%%%%%% extracting event
filtered_data = filtered_velocity;
event = event_func(filtered_data);
R_event = event;
R_Event_without_noise = [];
count = 1;
for i = 1:(length(R_event)-1)
    if R_event(i+1) - R_event(i) > 20
        R_Event_without_noise(count) = R_event(i);
        count = count + 1;
    end
end
%%%%%% plot
%plot(filtered_velocity)
%hold on 
%scatter(R_Event_without_noise, filtered_velocity(fix(R_Event_without_noise)))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Right Leg %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Trajectory extraction
LDM5_trajectory = Gait_sw_0001.Trajectories.Labeled.Data(find(strcmp(label, 'LDM5')),1, :);
LDM5_trajectory = reshape(LDM5_trajectory, [length(LDM5_trajectory), 1]);
%%%%%% velocity calculation
marker_trajectory = LDM5_trajectory;
velocity_result = velocity_func(marker_trajectory);
%%%%%% velocity normalized
velocity_result = velocity_result;
normalized_velocity = normalize(velocity_result);
%%%%%% velocity filtering
normalized_velocity = normalized_velocity;
filtered_velocity = filter_func(normalized_velocity);
%%%%%% extracting event
filtered_data = filtered_velocity;
event = event_func(filtered_data);
L_event = event;
L_Event_without_noise = [];
count = 1;
for i = 1:(length(L_event)-1)
    if L_event(i+1) - L_event(i) > 20
        L_Event_without_noise(count) = L_event(i);
        count = count + 1;
    end
end
%%%%%% plot
%plot(filtered_velocity)
%hold on 
%scatter(L_Event_without_noise, filtered_velocity(fix(L_Event_without_noise)))


