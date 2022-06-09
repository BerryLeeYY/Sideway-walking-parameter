%% load the data
clear
addpath 'C:\Users\a1003\OneDrive\桌面\Project_Review\sideward walking\inclinationa_angle\SW'
addpath 'C:\Users\a1003\OneDrive\桌面\Project_Review\sideward walking\inclinationa_angle\gait'
subject = ('S19_Gait_01_06');
subjfile = [(subject),'.mat'];
load(subjfile);

R_data = load(subject);
name =Gait_01_06;
label = Gait_01_06.Trajectories.Labeled.Labels;
%% analyze anterior-posterior inclination angle: coordination = 1
%  analyze medial-lateral incliantion angle: coordination = 2
coordination = 1 ;% 1: anterior-posterior, 2: medial-lateral, 3: up and down

%% store COP datat
FP1_data = name.Force(1).COP;
FP2_data = name.Force(2).COP;
FP3_data = name.Force(3).COP;
FP4_data = name.Force(4).COP;
FP5_data = name.Force(5).COP;
FP6_data = name.Force(6).COP;
FP7_data = name.Force(7).COP;


% store force data
FP1_data_Force = name.Force(1).Force;
FP2_data_Force = name.Force(2).Force;
FP3_data_Force = name.Force(3).Force;
FP4_data_Force = name.Force(4).Force;
FP5_data_Force = name.Force(5).Force;
FP6_data_Force = name.Force(6).Force;
FP7_data_Force = name.Force(7).Force;

data_len = length(FP1_data);
FP = zeros(6, data_len);

% corresponding time
frq = length(FP1_data) / length(name.Trajectories.Labeled.Data(26,1,:));
time = data_len / frq;

% form the COP matrix for 6 force plate

FP1_xyz = FP1_data(:,(1:frq:length(FP1_data)));
FP2_xyz = FP2_data(:,(1:frq:length(FP2_data)));
FP3_xyz = FP3_data(:,(1:frq:length(FP3_data)));
FP4_xyz = FP4_data(:,(1:frq:length(FP4_data)));
FP5_xyz = FP5_data(:,(1:frq:length(FP5_data)));
FP6_xyz = FP6_data(:,(1:frq:length(FP6_data)));
FP7_xyz = FP7_data(:,(1:frq:length(FP7_data)));



% form the force matrix for 6 force plate
FP1_xyz_Force = FP1_data_Force(:,(1:frq:length(FP1_data)));
FP2_xyz_Force = FP2_data_Force(:,(1:frq:length(FP2_data)));
FP3_xyz_Force = FP3_data_Force(:,(1:frq:length(FP3_data)));
FP4_xyz_Force = FP4_data_Force(:,(1:frq:length(FP4_data)));
FP5_xyz_Force = FP5_data_Force(:,(1:frq:length(FP5_data)));
FP6_xyz_Force = FP6_data_Force(:,(1:frq:length(FP6_data)));
FP7_xyz_Force = FP7_data_Force(:,(1:frq:length(FP7_data)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Right Leg %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Trajectory extraction
RDM5_trajectory = name.Trajectories.Labeled.Data(find(strcmp(label, 'RDM5')),1, :);
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
[event,  stride_event] = event_func(filtered_data);
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
LDM5_trajectory = name.Trajectories.Labeled.Data(find(strcmp(label, 'LDM5')),1, :);
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
[event,  stride_event] = event_func(filtered_data);
L_event = event;
L_Event_without_noise = [];
count = 1;
for i = 1:(length(L_event)-1)
    if L_event(i+1) - L_event(i) > 20
        L_Event_without_noise(count) = L_event(i);
        count = count + 1;
    end
end
%{
plot(filtered_velocity)
hold on 
scatter(L_Event_without_noise, filtered_velocity(fix(L_Event_without_noise)))
scatter(R_Event_without_noise,filtered_velocity(fix(R_Event_without_noise)))
%}


path = name.Trajectories.Labeled.Labels;
LCAL2_position = find(strcmp( path, 'LCAL2'));
RCAL2_position = find(strcmp( path, 'RCAL2'));

LCAL2_data = name.Trajectories.Labeled.Data(LCAL2_position,1:3,:);
RCAL2_data = name.Trajectories.Labeled.Data(RCAL2_position,1:3,:);
LCAL2_data_reshape = reshape(LCAL2_data, [3,time]);
RCAL2_data_reshape = reshape(RCAL2_data, [3,time]);

step_length = abs(RCAL2_data_reshape(1,:) - LCAL2_data_reshape(1,:));
step_width = abs(RCAL2_data_reshape(2,:) - LCAL2_data_reshape(2,:));
figure
plot(step_length)
hold on
scatter(L_Event_without_noise,step_length(fix(L_Event_without_noise)))
[peaks, loc] = findpeaks(step_length, "MinPeakProminence", 20);
findpeaks(step_length, "MinPeakProminence", 20);
loc


figure
plot(step_width)
hold on
scatter(L_Event_without_noise,step_width(fix(R_Event_without_noise)))

