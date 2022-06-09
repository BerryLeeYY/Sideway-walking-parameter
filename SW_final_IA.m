%% load the data
clear
addpath 'C:\Users\a1003\OneDrive\桌面\Project_Review\sideward walking\inclinationa_angle\SW'
addpath 'C:\Users\a1003\OneDrive\桌面\Project_Review\sideward walking\inclinationa_angle\gait'
addpath 'C:\Users\a1003\OneDrive\桌面\Project_Review\side_walk_new\SW'
subject = ('S09_Gait_sw_0001');
subjfile = [(subject),'.mat'];
load(subjfile);

R_data = load(subject);
name =Gait_sw_0001;
label = Gait_sw_0001.Trajectories.Labeled.Labels;
path = name.Trajectories.Labeled.Labels;
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
FP8_data = name.Force(8).COP;


% store force data
FP1_data_Force = name.Force(1).Force;
FP2_data_Force = name.Force(2).Force;
FP3_data_Force = name.Force(3).Force;
FP4_data_Force = name.Force(4).Force;
FP5_data_Force = name.Force(5).Force;
FP6_data_Force = name.Force(6).Force;
FP7_data_Force = name.Force(7).Force;
FP8_data_Force = name.Force(8).Force;

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
FP8_xyz = FP8_data(:,(1:frq:length(FP8_data)));



% form the force matrix for 6 force plate
FP1_xyz_Force = FP1_data_Force(:,(1:frq:length(FP1_data)));
FP2_xyz_Force = FP2_data_Force(:,(1:frq:length(FP2_data)));
FP3_xyz_Force = FP3_data_Force(:,(1:frq:length(FP3_data)));
FP4_xyz_Force = FP4_data_Force(:,(1:frq:length(FP4_data)));
FP5_xyz_Force = FP5_data_Force(:,(1:frq:length(FP5_data)));
FP6_xyz_Force = FP6_data_Force(:,(1:frq:length(FP6_data)));
FP7_xyz_Force = FP7_data_Force(:,(1:frq:length(FP7_data)));
FP8_xyz_Force = FP8_data_Force(:,(1:frq:length(FP8_data)));


%% Calculate the COM
%%%find the marker


LPSI_position = find(strcmp( path, 'LPSI'));
RPSI_position = find(strcmp( path, 'RPSI'));
LASI_position = find(strcmp( path, 'LASI'));
RASI_position = find(strcmp( path, 'RASI'));
LSHO_position = find(strcmp( path, 'LSHO'));
RSHO_position = find(strcmp( path, 'RSHO'));
LELL_position = find(strcmp( path, 'LELL'));
RELL_position = find(strcmp( path, 'RELL'));
LWRR_position = find(strcmp( path, 'LWRR'));
RWRR_position = find(strcmp( path, 'RWRR'));
LFLE_position = find(strcmp( path, 'LFLE'));
RFLE_position = find(strcmp( path, 'RFLE'));
LLMAL_position = find(strcmp( path, 'LLMAL'));
RLMAL_position = find(strcmp( path, 'RLMAL'));


%%
% 1: anterior-posterior, 2: medial-lateral, 3: up and down
LPSI_data = name.Trajectories.Labeled.Data(LPSI_position,1:3,:);
RPSI_data = name.Trajectories.Labeled.Data(RPSI_position,1:3,:);
LASI_data = name.Trajectories.Labeled.Data(LASI_position,1:3,:);
RASI_data = name.Trajectories.Labeled.Data(RASI_position,1:3,:);
LSHO_data = name.Trajectories.Labeled.Data(LSHO_position,1:3,:);
RSHO_data = name.Trajectories.Labeled.Data(RSHO_position,1:3,:);
LELL_data = name.Trajectories.Labeled.Data(LELL_position,1:3,:);
RELL_data = name.Trajectories.Labeled.Data(RELL_position,1:3,:);
%LWRR_data = name.Trajectories.Labeled.Data(LWRR_position,1:3,:);
%RWRR_data = name.Trajectories.Labeled.Data(RWRR_position,1:3,:);
LFLE_data = name.Trajectories.Labeled.Data(LFLE_position,1:3,:);
RFLE_data = name.Trajectories.Labeled.Data(RFLE_position,1:3,:);
LLMAL_data = name.Trajectories.Labeled.Data(LLMAL_position,1:3,:);
RLMAL_data = name.Trajectories.Labeled.Data(RLMAL_position,1:3,:);

%%
LPSI_data = reshape(LPSI_data, [3, time]);
RPSI_data = reshape(RPSI_data, [3,time]);
LASI_data = reshape(LASI_data, [3,time]);
RASI_data = reshape(RASI_data, [3,time]);
LSHO_data = reshape(LSHO_data, [3,time]);
RSHO_data = reshape(RSHO_data, [3,time]);
LELL_data = reshape(LELL_data, [3,time]);
RELL_data = reshape(RELL_data, [3,time]);
%LWRR_data = reshape(LWRR_data, [3,time]);
%RWRR_data = reshape(RWRR_data, [3,time]);
LFLE_data = reshape(LFLE_data, [3,time]);
RFLE_data = reshape(RFLE_data, [3,time]);
LLMAL_data = reshape(LLMAL_data, [3,time]);
RLMAL_data = reshape(RLMAL_data, [3,time]);


%%
LASI = LASI_data;
LPSI = LPSI_data;
RASI = RASI_data;
RPSI = RPSI_data;

% calculate the hip marker position
[hip_center, L_hip_center, R_hip_center] = hip_markers(LASI, LPSI, RASI, RPSI);

% store the position data from each marker
L_shoulder  = LSHO_data;
R_shoulder  = RSHO_data;
L_elbow     = LELL_data;
R_elbow     = RELL_data;
L_hand      = "missing_marker";
R_hand	    = "missing_marker";
L_knee	    = LFLE_data;
R_knee      = RFLE_data;
L_ankle	    = LLMAL_data;
R_ankle     = RLMAL_data;
hip_center = hip_center;
L_hip_center = L_hip_center;
R_hip_center = R_hip_center;


New_COM = COM_function(time, L_shoulder, R_shoulder, L_elbow, R_elbow, L_hand, R_hand, L_knee, R_knee, L_ankle, R_ankle,hip_center, L_hip_center, R_hip_center, 1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Right Leg %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Trajectory extraction
RDM5_trajectory = name.Trajectories.Labeled.Data(find(strcmp(label, 'RDM5')),1, :);
RDM5_trajectory = reshape(RDM5_trajectory, [length(RDM5_trajectory), 1]);
%%%%%% velocity calculation
marker_trajectory = RDM5_trajectory;
velocity_result = velocity_func(marker_trajectory);
%%%%%% velocity normalized
velocity_result = velocity_result(:,1);
filtered_velocity = SW_filter_func(velocity_result);
%%%%%% extracting event
filtered_data = filtered_velocity;
[event,  stride_event] = event_func(filtered_data);
R_event = event;
R_Event_without_noise = [];
count = 1;
R_foot_contact_event = R_event(stride_event == 1);
R_foot_off_event = R_event(stride_event == 0);
R_foot_off_event = R_foot_off_event(2:end);
for i = 1:(length(R_event(stride_event == 1)))
    try
        if (R_foot_contact_event(i+1) - R_foot_contact_event(i)) > 20
            R_Event_without_noise(count) = R_foot_contact_event(i);
            count = count + 1;
        end
    catch
        continue
    end
end
%%%%%% plot
figure
plot(filtered_velocity)
hold on 
scatter(R_foot_contact_event, filtered_velocity(fix(R_foot_contact_event)))
title("Right foot contact")

figure
plot(filtered_velocity)
hold on 
scatter(R_foot_off_event, filtered_velocity(fix(R_foot_off_event)))
title("Right foot off")



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Left Leg %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Trajectory extraction
LDM5_trajectory = name.Trajectories.Labeled.Data(find(strcmp(label, 'LDM5')),1, :);
LDM5_trajectory = reshape(LDM5_trajectory, [length(LDM5_trajectory), 1]);
%%%%%% velocity calculation
marker_trajectory = LDM5_trajectory;
velocity_result = velocity_func(marker_trajectory);
%%%%%% velocity normalized
velocity_result = velocity_result(:,1);
filtered_velocity = SW_filter_func(velocity_result);
%%%%%% extracting event
filtered_data = filtered_velocity;
[event,  stride_event] = event_func(filtered_data);
L_event = event;
L_Event_without_noise = [];
count = 1;
L_foot_contact_event = L_event(stride_event == 1);
L_foot_off_event = L_event(stride_event == 0);
L_foot_off_event = L_foot_off_event(2:end);
for i = 1:(length(L_event(stride_event == 1)))
    try
        if (L_foot_contact_event(i+1) - L_foot_contact_event(i)) > 20
            L_Event_without_noise(count) = L_foot_contact_event(i);
            count = count + 1;
        end
    catch
        continue
    end
end

%%%%%% plot
figure
plot(filtered_velocity)
hold on 
scatter(L_foot_contact_event, filtered_velocity(fix(L_foot_contact_event)))
title("Left foot contact")

figure
plot(filtered_velocity)
hold on 
scatter(L_foot_off_event, filtered_velocity(fix(L_foot_off_event)))
title("Light foot off")


%%%%%%%%%%%%%%%  filter stepping force platform  %%%%%%%%%%%%%%%
force_plat_data = zeros(7, length(FP4_xyz(1,:)));
force_plat_data(1,:) = FP1_xyz(coordination,:);
force_plat_data(2,:) = FP2_xyz(coordination,:);
force_plat_data(3,:) = FP3_xyz(coordination,:);
force_plat_data(4,:) = FP4_xyz(coordination,:);
force_plat_data(5,:) = FP5_xyz(coordination,:);
force_plat_data(6,:) = FP6_xyz(coordination,:);
force_plat_data(7,:) = FP7_xyz(coordination,:);
force_plat_data(8,:) = FP8_xyz(coordination,:);

force_plat_forcedata = zeros(7, length(FP4_xyz(1,:)));
force_plat_forcedata(1,:) = FP1_xyz_Force(3,:);
force_plat_forcedata(2,:) = FP2_xyz_Force(3,:);
force_plat_forcedata(3,:) = FP3_xyz_Force(3,:);
force_plat_forcedata(4,:) = FP4_xyz_Force(3,:);
force_plat_forcedata(5,:) = FP5_xyz_Force(3,:);
force_plat_forcedata(6,:) = FP6_xyz_Force(3,:);
force_plat_forcedata(7,:) = FP7_xyz_Force(3,:);
force_plat_forcedata(8,:) = FP8_xyz_Force(3,:);


force_plat_identification = zeros(7, length(FP4_xyz(1,:)));
force_plat_identification(1,:) = FP1_xyz(1,:);
force_plat_identification(2,:) = FP2_xyz(1,:);
force_plat_identification(3,:) = FP3_xyz(1,:);
force_plat_identification(4,:) = FP4_xyz(1,:);
force_plat_identification(5,:) = FP5_xyz(1,:);
force_plat_identification(6,:) = FP6_xyz(1,:);
force_plat_identification(7,:) = FP7_xyz(1,:);
force_plat_identification(8,:) = FP8_xyz(1,:);

force_plat_valid_signal = zeros(1,length(FP4_xyz(1,:)));
force_data = zeros(1,length(FP4_xyz(1,:)));
count = 1;
for i = 1:8
    if std(force_plat_forcedata(i,:)) > 100
        force_plat_valid_signal(count, :) = force_plat_data(i,:);
        force_data(count, :) = force_plat_forcedata(i, :);
        count = count + 1;
    end
end
disp(size(force_plat_valid_signal))
disp(size(force_data))
%% COP Processing


for i = 1:size(force_plat_valid_signal,1)
    
    for ii = 1:size(force_plat_valid_signal,2)
        if abs(force_data(i,ii)) > 10 %(10~50)
            x = force_plat_valid_signal(i,ii);     
        else
            x = 0;
        end
        y(i, ii) = x;
    end
end

if  size(force_plat_valid_signal,1) == 3
    figure
    plot(y(1,:))
    hold on
    plot(y(2,:))
    plot(y(3,:))
    title('COP trajectory in each force platform')
elseif size(force_plat_valid_signal,1) == 2
    figure
    plot(y(1,:))
    hold on
    plot(y(2,:))
    title('COP trajectory in each force platform')
elseif size(force_plat_valid_signal,1) == 4
    figure
    plot(y(1,:))
    hold on
    plot(y(2,:))
    plot(y(3,:))
    plot(y(4,:))
    title('COP trajectory in each force platform')
elseif size(force_plat_valid_signal,1) == 5
    figure
    plot(y(1,:))
    hold on
    plot(y(2,:))
    plot(y(3,:))
    plot(y(4,:))
    plot(y(5,:))
    title('COP trajectory in each force platform')
elseif size(force_plat_valid_signal,1) == 6
    figure
    plot(y(1,:))
    hold on
    plot(y(2,:))
    plot(y(3,:))
    plot(y(4,:))
    plot(y(5,:))
    plot(y(6,:))
    title('COP trajectory in each force platform')
elseif size(force_plat_valid_signal,1) == 7
    figure
    plot(y(1,:))
    hold on
    plot(y(2,:))
    plot(y(3,:))
    plot(y(4,:))
    plot(y(5,:))
    plot(y(6,:))
    plot(y(7,:))
    title('COP trajectory in each force platform')
end

if size(force_plat_valid_signal, 1) == 3
    for i = 1:size(y,2)
        if y(1,i) ~= 0 && y(2,i) ~= 0
            value = (y(1,i)*force_data(1,i) + y(2,i)*force_data(2,i)) / (force_data(1,i) + force_data(2,i));
        elseif y(1,i) ~= 0 && y(3,i) ~= 0
            value = (y(1,i)*force_data(1,i) + y(3,i)*force_data(3,i)) / (force_data(1,i) + force_data(3,i));
        elseif y(2,i) ~= 0 && y(3,i) ~= 0
            value = (y(2,i)*force_data(2,i) + y(3,i)*force_data(3,i)) / (force_data(2,i) + force_data(3,i));
        elseif y(1,i) ~= 0 && y(2,i) == 0 && y(3,i) == 0
            value = y(1,i);
        elseif y(1,i) == 0 && y(2,i) ~= 0 && y(3,i) == 0
            value = y(2,i);
        elseif y(1,i) == 0 && y(2,i) == 0 && y(3,i) ~= 0
            value = y(3,i);
        else
            value = 0;
        end
        final_cop(i) = value;
    end

elseif size(force_plat_valid_signal, 1) == 2
    for i = 1:size(y,2)
        if y(1,i) ~= 0 && y(2,i) ~= 0
            value = (y(1,i)*force_data(1,i) + y(2,i)*force_data(2,i)) / (force_data(1,i) + force_data(2,i));
        elseif y(1,i) ~= 0 && y(2,i) == 0 
            value = y(1,i);
        elseif y(1,i) == 0 && y(2,i) ~= 0 
            value = y(2,i);
        else
            value = 0;
        end
        final_cop(i) = value;

    end 
    

elseif size(force_plat_valid_signal, 1) == 4
    for i = 1:size(y,2)
        if y(1,i) ~= 0 && y(2,i) ~= 0
            value = (y(1,i)*force_data(1,i) + y(2,i)*force_data(2,i)) / (force_data(1,i) + force_data(2,i));
        elseif y(1,i) ~= 0 && y(3,i) ~= 0
            value = (y(1,i)*force_data(1,i) + y(3,i)*force_data(3,i)) / (force_data(1,i) + force_data(3,i));
        elseif y(2,i) ~= 0 && y(3,i) ~= 0
            value = (y(2,i)*force_data(2,i) + y(3,i)*force_data(3,i)) / (force_data(2,i) + force_data(3,i));
        elseif y(1,i) ~= 0 && y(4,i) ~= 0
            value = (y(1,i)*force_data(1,i) + y(4,i)*force_data(4,i)) / (force_data(1,i) + force_data(4,i));
        elseif y(2,i) ~= 0 && y(4,i) ~= 0
            value = (y(2,i)*force_data(2,i) + y(4,i)*force_data(4,i)) / (force_data(2,i) + force_data(4,i));
        elseif y(3,i) ~= 0 && y(4,i) ~= 0
            value = (y(3,i)*force_data(3,i) + y(4,i)*force_data(4,i)) / (force_data(3,i) + force_data(4,i));    
            
        elseif y(1,i) ~= 0 && y(2,i) == 0 && y(3,i) == 0 && y(4,i) == 0
            value = y(1,i);
        elseif y(1,i) == 0 && y(2,i) ~= 0 && y(3,i) == 0 && y(4,i) == 0
            value = y(2,i);
        elseif y(1,i) == 0 && y(2,i) == 0 && y(3,i) ~= 0 && y(4,i) == 0
            value = y(3,i);
        elseif y(1,i) == 0 && y(2,i) == 0 && y(3,i) == 0 && y(4,i) ~= 0
            value = y(4,i);
        else
            value = 0;
        end
        final_cop(i) = value;
    end
    
elseif size(force_plat_valid_signal, 1) == 5
    for i = 1:size(y,2)
        if y(1,i) ~= 0 && y(2,i) ~= 0
            value = (y(1,i)*force_data(1,i) + y(2,i)*force_data(2,i)) / (force_data(1,i) + force_data(2,i));
        elseif y(1,i) ~= 0 && y(3,i) ~= 0
            value = (y(1,i)*force_data(1,i) + y(3,i)*force_data(3,i)) / (force_data(1,i) + force_data(3,i));
        elseif y(1,i) ~= 0 && y(4,i) ~= 0
            value = (y(1,i)*force_data(1,i) + y(4,i)*force_data(4,i)) / (force_data(1,i) + force_data(4,i));
        elseif y(1,i) ~= 0 && y(5,i) ~= 0
            value = (y(1,i)*force_data(1,i) + y(5,i)*force_data(5,i)) / (force_data(1,i) + force_data(5,i));
            
        elseif y(2,i) ~= 0 && y(3,i) ~= 0
            value = (y(2,i)*force_data(2,i) + y(3,i)*force_data(3,i)) / (force_data(2,i) + force_data(3,i));
        elseif y(2,i) ~= 0 && y(4,i) ~= 0
            value = (y(2,i)*force_data(2,i) + y(4,i)*force_data(4,i)) / (force_data(2,i) + force_data(4,i));
        elseif y(2,i) ~= 0 && y(5,i) ~= 0
            value = (y(2,i)*force_data(2,i) + y(5,i)*force_data(4,i)) / (force_data(2,i) + force_data(5,i));
            
        elseif y(3,i) ~= 0 && y(4,i) ~= 0
            value = (y(3,i)*force_data(3,i) + y(4,i)*force_data(4,i)) / (force_data(3,i) + force_data(4,i));
        elseif y(3,i) ~= 0 && y(5,i) ~= 0
            value = (y(3,i)*force_data(3,i) + y(5,i)*force_data(5,i)) / (force_data(3,i) + force_data(5,i)); 
        
        elseif y(4,i) ~= 0 && y(5,i) ~= 0
            value = (y(4,i)*force_data(4,i) + y(5,i)*force_data(5,i)) / (force_data(4,i) + force_data(5,i)); 
            
        elseif y(1,i) ~= 0 && y(2,i) == 0 && y(3,i) == 0 && y(4,i) == 0
            value = y(1,i);
        elseif y(1,i) == 0 && y(2,i) ~= 0 && y(3,i) == 0 && y(4,i) == 0
            value = y(2,i);
        elseif y(1,i) == 0 && y(2,i) == 0 && y(3,i) ~= 0 && y(4,i) == 0
            value = y(3,i);
        elseif y(1,i) == 0 && y(2,i) == 0 && y(3,i) == 0 && y(4,i) ~= 0
            value = y(4,i);
            
        elseif y(1,i) ~= 0 && y(2,i) == 0 && y(3,i) == 0 && y(4,i) == 0 && y(5,i) == 0
            value = y(1,i);
        elseif y(1,i) == 0 && y(2,i) ~= 0 && y(3,i) == 0 && y(4,i) == 0 && y(5,i) == 0
            value = y(2,i);
        elseif y(1,i) == 0 && y(2,i) == 0 && y(3,i) ~= 0 && y(4,i) == 0 && y(5,i) == 0
            value = y(3,i);
        elseif y(1,i) == 0 && y(2,i) == 0 && y(3,i) == 0 && y(4,i) ~= 0 && y(5,i) == 0
            value = y(4,i);
        elseif y(1,i) == 0 && y(2,i) == 0 && y(3,i) == 0 && y(4,i) ~= 0 && y(5,i) == 0
            value = y(5,i);
        
       
        else
            value = 0;
        end
        final_cop(i) = value;
    end
elseif size(force_plat_valid_signal, 1) == 6
    for i = 1:size(y,2)
        if y(1,i) ~= 0 && y(2,i) ~= 0
            value = (y(1,i)*force_data(1,i) + y(2,i)*force_data(2,i)) / (force_data(1,i) + force_data(2,i));
        elseif y(1,i) ~= 0 && y(3,i) ~= 0
            value = (y(1,i)*force_data(1,i) + y(3,i)*force_data(3,i)) / (force_data(1,i) + force_data(3,i));
        elseif y(1,i) ~= 0 && y(4,i) ~= 0
            value = (y(1,i)*force_data(1,i) + y(4,i)*force_data(4,i)) / (force_data(1,i) + force_data(4,i));
        elseif y(1,i) ~= 0 && y(5,i) ~= 0
            value = (y(1,i)*force_data(1,i) + y(5,i)*force_data(5,i)) / (force_data(1,i) + force_data(5,i));
        elseif y(1,i) ~= 0 && y(6,i) ~= 0
            value = (y(1,i)*force_data(1,i) + y(6,i)*force_data(6,i)) / (force_data(1,i) + force_data(6,i));
            
        elseif y(2,i) ~= 0 && y(3,i) ~= 0
            value = (y(2,i)*force_data(2,i) + y(3,i)*force_data(3,i)) / (force_data(2,i) + force_data(3,i));
        elseif y(2,i) ~= 0 && y(4,i) ~= 0
            value = (y(2,i)*force_data(2,i) + y(4,i)*force_data(4,i)) / (force_data(2,i) + force_data(4,i));
        elseif y(2,i) ~= 0 && y(5,i) ~= 0
            value = (y(2,i)*force_data(2,i) + y(5,i)*force_data(5,i)) / (force_data(2,i) + force_data(5,i));
        elseif y(2,i) ~= 0 && y(6,i) ~= 0
            value = (y(2,i)*force_data(2,i) + y(6,i)*force_data(6,i)) / (force_data(2,i) + force_data(6,i));
            
        elseif y(3,i) ~= 0 && y(4,i) ~= 0
            value = (y(3,i)*force_data(3,i) + y(4,i)*force_data(4,i)) / (force_data(3,i) + force_data(4,i));
        elseif y(3,i) ~= 0 && y(5,i) ~= 0
            value = (y(3,i)*force_data(3,i) + y(5,i)*force_data(5,i)) / (force_data(3,i) + force_data(5,i)); 
        elseif y(3,i) ~= 0 && y(6,i) ~= 0
            value = (y(3,i)*force_data(3,i) + y(6,i)*force_data(6,i)) / (force_data(3,i) + force_data(6,i)); 
        
        elseif y(4,i) ~= 0 && y(5,i) ~= 0
            value = (y(4,i)*force_data(4,i) + y(5,i)*force_data(5,i)) / (force_data(4,i) + force_data(5,i)); 
        elseif y(4,i) ~= 0 && y(6,i) ~= 0
            value = (y(4,i)*force_data(4,i) + y(6,i)*force_data(6,i)) / (force_data(4,i) + force_data(6,i)); 
            
        elseif y(5,i) ~= 0 && y(6,i) ~= 0
            value = (y(5,i)*force_data(5,i) + y(6,i)*force_data(6,i)) / (force_data(5,i) + force_data(6,i)); 
     
            
        elseif y(1,i) ~= 0 && y(2,i) == 0 && y(3,i) == 0 && y(4,i) == 0 && y(5,i) == 0 && y(6,i) == 0
            value = y(1,i);
        elseif y(1,i) == 0 && y(2,i) ~= 0 && y(3,i) == 0 && y(4,i) == 0 && y(5,i) == 0 && y(6,i) == 0
            value = y(2,i);
        elseif y(1,i) == 0 && y(2,i) == 0 && y(3,i) ~= 0 && y(4,i) == 0 && y(5,i) == 0 && y(6,i) == 0
            value = y(3,i);
        elseif y(1,i) == 0 && y(2,i) == 0 && y(3,i) == 0 && y(4,i) ~= 0 && y(5,i) == 0 && y(6,i) == 0
            value = y(4,i);
        elseif y(1,i) == 0 && y(2,i) == 0 && y(3,i) == 0 && y(4,i) == 0 && y(5,i) ~= 0 && y(6,i) == 0
            value = y(5,i);
        elseif y(1,i) == 0 && y(2,i) == 0 && y(3,i) == 0 && y(4,i) == 0 && y(5,i) == 0 && y(6,i) ~= 0
            value = y(6,i);
        
       
        else
            value = 0;
        end
        final_cop(i) = value;
    end
  
    
end
    figure
    plot(final_cop)
    title('Combination of COP trajectory')

   
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  IA calculation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COM = New_COM;
threshold = 30;
IA = nan(size(force_plat_valid_signal, 1), length(force_plat_valid_signal));

for i = 1:length(final_cop)
    if abs(final_cop(i)) > 0
        ia(i) =  (final_cop(i)-New_COM(coordination,i))/New_COM(3,i);
    else
        ia(i) =  nan;
    end
end
IA = atand(ia);

%{
for i  = 1:size(force_plat_valid_signal, 1)
    IA(i,:) = filter_IA(force_plat_valid_signal(i,:), COM, threshold);
end
IA_all = zeros(1, length(IA));
for i = 1:size(IA, 1)
    x = IA(i, :);
    IA_all = IA_all + x;
end
%IA_all = IA_1 + IA_2 + IA_3;
%}
figure
plot(IA)
hold on
scatter(L_foot_contact_event, IA(fix(L_foot_contact_event)));
scatter(R_foot_contact_event, IA(fix(R_foot_contact_event)));
title('Inclination angle with gait event (foot contact)')

figure
plot(IA)
hold on
scatter(L_foot_off_event, IA(fix(L_foot_off_event)));
scatter(R_foot_off_event, IA(fix(R_foot_off_event)));
title('Inclination angle with gait event (foot contact)')



%% step length and step size
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
scatter(R_foot_contact_event, step_length(fix(R_foot_contact_event)));
title('Step length with gait event')

FS_RE_IA = IA(fix(R_foot_contact_event(step_length(fix(R_foot_contact_event)) > 300)));
FS_LE_IA = IA(fix(L_foot_contact_event(step_length(fix(L_foot_contact_event)) > 300)));
FS_R_e = R_foot_contact_event(step_length(fix(R_foot_contact_event)) > 300);
FS_L_e = L_foot_contact_event(step_length(fix(L_foot_contact_event)) > 300);
figure
plot(IA)
hold on
scatter(FS_L_e, IA(fix(FS_L_e)));
scatter(FS_R_e, IA(fix(FS_R_e)));
title('Inclination angle with gait event (foot contact widest)')


FO_RE_IA = IA(fix(R_foot_off_event(step_length(fix(R_foot_off_event)) > 300)));
FO_LE_IA = IA(fix(L_foot_off_event(step_length(fix(L_foot_off_event)) > 300)));
FO_R_e = R_foot_off_event(step_length(fix(R_foot_off_event)) > 300);
FO_L_e = L_foot_off_event(step_length(fix(L_foot_off_event)) > 300);

figure
plot(IA)
hold on
scatter(FO_L_e, IA(fix(FO_L_e)));
scatter(FO_R_e, IA(fix(FO_R_e)));
title('Inclination angle with gait event (foot off widest)')
%{
RE_SL = step_length(fix(R_e));
LE_SL = step_length(fix(L_e));
count = 1;
for i = 1:length(LE_SL)
    if 1
        L_Peak_SL(count) = abs(LE_SL(1,i));
        count = count + 1;
    end
end
count = 1;
for i = 1:length(RE_SL)
    if 1
        R_Peak_SL(count) = abs(RE_SL(1,i));
        count = count + 1;
    end
end
        
Peak_SL = [L_Peak_SL, R_Peak_SL];

figure
plot(step_width)
hold on
%scatter(L_Event_without_noise,step_width(fix(R_Event_without_noise)))
title('Step width with gait event')

RE_SW = step_width(fix(R_e));
LE_SW = step_width(fix(L_e));
count = 1;
for i = 1:length(LE_SW)
    if 1
        L_Peak_SW(count) = abs(LE_SW(1,i));
        count = count + 1;
    end
end
count = 1;
for i = 1:length(RE_SW)
    if 1
        R_Peak_SW(count) = abs(RE_SW(1,i));
        count = count + 1;
    end
end
        
Peak_SW = [L_Peak_SW, R_Peak_SW];

COM_x = COM(1,:);
COM_y = COM(2,:);
figure
corrcoef(COM_x(fix(L_Event_without_noise)),IA(fix(L_Event_without_noise)), 'rows','complete')
scatter(COM_x, IA)
scatter(step_length, IA)

RE_COM_x = COM_x(fix(R_e));
LE_COM_x = COM_x(fix(L_e));
count = 1;
for i = 1:length(LE_COM_x)
    if 1
        L_Peak_COM_x(count) = abs(LE_COM_x(1,i));
        count = count + 1;
    end
end
count = 1;
for i = 1:length(RE_COM_x)
    if 1
        R_Peak_COM_x(count) = abs(RE_COM_x(1,i));
        count = count + 1;
    end
end
        
Peak_COM_x = [L_Peak_COM_x, R_Peak_COM_x];

RE_COM_y = COM_y(fix(R_e));
LE_COM_y = COM_y(fix(L_e));
count = 1;
for i = 1:length(LE_COM_y)
    if 1
        L_Peak_COM_y(count) = abs(LE_COM_y(1,i));
        count = count + 1;
    end
end
count = 1;
for i = 1:length(RE_COM_x)
    if 1
        R_Peak_COM_y(count) = abs(RE_COM_y(1,i));
        count = count + 1;
    end
end

Peak_COM_y = [L_Peak_COM_y, R_Peak_COM_y];
%}
%{
if coordination == 1
    corr = corrcoef(Peak_IA, Peak_SL, 'rows','complete');
    disp("The correlation between IA and step length is " + num2str(corr(1,2)))
    figure
    scatter(Peak_IA, Peak_SL)
    title('Peak IA vs Peak step_length')
    corr_COM = corrcoef(Peak_IA, Peak_COM_x, 'rows','complete');
    disp("The correlation between IA and COM is " + num2str(corr_COM(1,2)))
    figure
    scatter(Peak_IA, Peak_COM_x)
    title('Peak IA frontal vs Peak COM position')
    csvwrite([subject(1:4) + "peak_frontal_IA.csv"], Peak_IA)
    csvwrite([subject(1:4) + "peak_SL.csv"], Peak_SL)
    csvwrite([subject(1:4) + "peak_SW.csv"], Peak_SW)
    csvwrite([subject(1:4) + "peak_COM_x.csv"], Peak_COM_x)
    csvwrite([subject(1:4) + "peak_COM_y.csv"], Peak_COM_y)
    csvwrite([subject(1:4) + "peak_frontal_all.csv"], [Peak_IA; Peak_SL; Peak_SW; Peak_COM_x; Peak_COM_y])
    
elseif coordination == 2
    corr = corrcoef(Peak_IA, Peak_SW, 'rows','complete');
    disp("The correlation between IA and step width is "  + num2str(corr(1,2)))
    figure
    scatter(Peak_IA, Peak_SW)
    title('Peak IA vs Peak step width')
    corr_COM = corrcoef(Peak_IA, Peak_COM_y, 'rows','complete');
    disp("The correlation between IA and COM is " + num2str(corr_COM(1,2)))
    figure
    scatter(Peak_IA, Peak_COM_y)
    title('Peak IA sagittal vs Peak COM position')
    csvwrite([subject(1:4) + "peak_sagittal_IA.csv"], Peak_IA)
    csvwrite([subject(1:4) + "peak_SL.csv"], Peak_SL)
    csvwrite([subject(1:4) + "peak_SW.csv"], Peak_SW)
    csvwrite([subject(1:4) + "peak_COM_x.csv"], Peak_COM_x)
    csvwrite([subject(1:4) + "peak_COM_y.csv"], Peak_COM_y)
    csvwrite([subject(1:4) + "peak_sagittal_all.csv"], [Peak_IA; Peak_SL; Peak_SW; Peak_COM_x; Peak_COM_y])
end
%}


