%%
%%% refer article: Ohtsu 2019 Investigation of balance strategy over gait cycle based on margin of stability
%%% refer article: https://royalsocietypublishing.org/doi/pdf/10.1098/rsif.2020.0194
%%% velocity -> using the window to avoid the jumping 
%%% position unit
%%% check through the video
%% load the data
clear
addpath 'C:\Users\a1003\OneDrive\桌面\Project_Review\sideward walking\inclinationa_angle\gait'
addpath 'C:\Users\a1003\OneDrive\桌面\Project_Review\side_walk_new\SW'
subject = ('S01_Gait_0001_2');
subjfile = [(subject),'.mat'];
load(subjfile);
%data_label = 12;
R_data = load(subject);
name =Gait_0001___2;
label = Gait_0001___2.Trajectories.Labeled.Labels;

FP1_data = name.Force(1).COP;
data_len = length(FP1_data);
FP = zeros(6, data_len);

% corresponding time
frq = length(FP1_data) / length(name.Trajectories.Labeled.Data(26,1,:));
time = data_len / frq;
%% COM calculation
%%% Calculate the COM
%%%find the marker
path = name.Trajectories.Labeled.Labels;

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


%%%
%coordination = 2 ;% 1: anterior-posterior, 2: medial-lateral, 3: up and down
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


%%%
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


%%%
LASI = LASI_data;
LPSI = LPSI_data;
RASI = RASI_data;
RPSI = RPSI_data;

[hip_center, L_hip_center, R_hip_center] = hip_markers(LASI, LPSI, RASI, RPSI);

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
%% XCOM
% w0 = sqrt(g / l)
% xCOM = COM + Vcom / w0
% l = pendulum length, leg length
COM = New_COM ./ 1000;
dt = 1 / 200;


%%% the code below is to form the velocity data in different axis for different marker
% SHO = axis_v_data(1,:), ELL = axis_v_data(2,:), ELM = axis_v_data(3,:), WRR = axis_v_data(4,:), WRU = axis_v_data(5,:)


x_v_data = zeros(1, time );
for i = 1:(time-1)
    x_v_data(1,i) = ((COM(1, i+1) - COM(1, i)) / dt);
end


y_v_data = zeros(1, time );
for i = 1:(time-1)
    y_v_data(1,i) = ((COM(2, i+1) - COM(2, i)) / dt);
end

z_v_data = zeros(1, time );
for i = 1:(time-1)
    z_v_data(1,i) = ((COM(3, i+1) - COM(3, i)) / dt);
end

total = zeros(3,time);
for i = 1:(time-1)
    total(1,i) = x_v_data(1,i);
    total(2,i) = y_v_data(1,i);
    total(3,i) = z_v_data(1,i);
end

% plot the velocity
%plot(total(1,:))
%hold on 
% plot the trajectory
%plot(WRR(3,:))


duration = 10;
g = 9.81;
l = (nanmean(COM(3,:)));
w0 = sqrt( g / l);
Vcom = x_v_data(1,:);
XCOM = (COM) + (total/ w0);
%% Marker information (position and velocity
RCAL1_position = find(strcmp( path, 'RCAL1'));
RCAL1_data = name.Trajectories.Labeled.Data(RCAL1_position,1,:);
RCAL1_data = reshape(RCAL1_data, [time, 1]);
RCAL1_height = name.Trajectories.Labeled.Data(RCAL1_position,3,:);
RCAL1_height = reshape(RCAL1_height, [1, time]);

LCAL1_position = find(strcmp( path, 'LCAL1'));
LCAL1_data = name.Trajectories.Labeled.Data(LCAL1_position,1,:);
LCAL1_data = reshape(LCAL1_data, [time, 1]);
LCAL1_height = name.Trajectories.Labeled.Data(LCAL1_position,3,:);
LCAL1_height = reshape(LCAL1_height, [1, time]);

RDM2_position = find(strcmp( path, 'RDM2'));
RDM2_data = name.Trajectories.Labeled.Data(RDM2_position,1,:);
RDM2_data = reshape(RDM2_data, [time, 1]);
RDM2_height = name.Trajectories.Labeled.Data(RDM2_position,3,:);
RDM2_height = reshape(RDM2_height, [1, time]);

LDM2_position = find(strcmp( path, 'LDM2'));
LDM2_data = name.Trajectories.Labeled.Data(LDM2_position,1,:);
LDM2_data = reshape(LDM2_data, [time, 1]);
LDM2_height = name.Trajectories.Labeled.Data(LDM2_position,3,:);
LDM2_height = reshape(LDM2_height, [1, time]);

marker_trajectory = RCAL1_data;
RCAL1_data_velocity_result = abs(velocity_func(marker_trajectory))/10;
marker_trajectory = LCAL1_data;
LCAL1_data_velocity_result = abs(velocity_func(marker_trajectory))/10;
marker_trajectory = RDM2_data;
RDM2_data_velocity_result = abs(velocity_func(marker_trajectory))/10;
marker_trajectory = LDM2_data;
LDM2_data_velocity_result = abs(velocity_func(marker_trajectory))/10;

%% Identification of the event
velocity_threshold = 100;
CAL1_position_threshold = 30;
DM2_position_threshold = 40;
%%% Heel Strike to before foot flat (HS_BFF) %%%
%%% velocity of L(R)CAL1 below velocity_threshold  
%%% velocity of R(L)CAL1 above velocity_threshold  
%%% velocity of R(L)DM2 below velocity_threshold   
%%% velocity of L(R)DM2 above velocity_threshold   
%%% L(R)CAL1 below CAL1_position_threshold                 
%%% R(L)CAL1 above CAL1_position_threshold                 
%%% R(L)DM2 below DM2_position_threshold                  
%%% L(R)DM2 above DM2_position_threshold                  

%Velocity rules
%RCAL1_data_velocity_result(i,1) < velocity_threshold && LCAL1_data_velocity_result(i,1) > velocity_threshold && RDM2_data_velocity_result(i,1) < velocity_threshold && LDM2_data_velocity_result(i,1) > velocity_threshold 
%LCAL1_data_velocity_result(i,1) < velocity_threshold && RCAL1_data_velocity_result(i,1) > velocity_threshold && LDM2_data_velocity_result(i,1) < velocity_threshold && RDM2_data_velocity_result(i,1) > velocity_threshold

%Position rules
%LCAL1_height(1, i) < CAL1_position_threshold && RCAL1_height(1, i) > CAL1_position_threshold && RDM2_height(1, i) < DM2_position_threshold && LDM2_height(1, i) > DM2_position_threshold
%RCAL1_height(1, i) < CAL1_position_threshold && LCAL1_height(1, i) > CAL1_position_threshold && LDM2_height(1, i) < DM2_position_threshold && RDM2_height(1, i) > DM2_position_threshold


for i = 1:time
    if LCAL1_height(1, i) < CAL1_position_threshold && RCAL1_height(1, i) > CAL1_position_threshold && RDM2_height(1, i) < DM2_position_threshold && LDM2_height(1, i) > DM2_position_threshold
        HS_BFF(i) = 1;
    elseif  RCAL1_height(1, i) < CAL1_position_threshold && LCAL1_height(1, i) > CAL1_position_threshold && LDM2_height(1, i) < DM2_position_threshold && RDM2_height(1, i) > DM2_position_threshold
        HS_BFF(i) = 1;
    else 
        HS_BFF(i) = 0;
    end
end

%%% Foot flat to before heel off (FF_BHO) %%%
%%% velocity of L(R)CAL1 above 150 mm/s  LCAL1_data_velocity_result(i,1) > 150
%%% velocity of R(L)CAL1 below 150 mm/s  RCAL1_data_velocity_result(i,1) < 150
%%% velocity of R(L)DM2 below 150 mm/s   RDM2_data_velocity_result(i,1) < 150 
%%% velocity of L(R)DM2 above 150 mm/s   LDM2_data_velocity_result(i,1) > 150 
%%% L(R)CAL1 above 35 mm                 LCAL1_height(1, i) > 35 
%%% R(L)CAL1 below 35 mm                 RCAL1_height(1, i) < 35
%%% R(L)DM2 below 35 mm                  RDM2_height(1, i) < 35
%%% L(R)DM2 above 35 mm                  LDM2_height(1, i) > 35

%Velocity rules
%LCAL1_data_velocity_result(i,1) > velocity_threshold && RCAL1_data_velocity_result(i,1) < velocity_threshold && RDM2_data_velocity_result(i,1) < velocity_threshold && LDM2_data_velocity_result(i,1) > velocity_threshold
%RCAL1_data_velocity_result(i,1) > velocity_threshold && LCAL1_data_velocity_result(i,1) < velocity_threshold && LDM2_data_velocity_result(i,1) < velocity_threshold && RDM2_data_velocity_result(i,1) > velocity_threshold

%Position rules
%LCAL1_height(1, i) > CAL1_position_threshold && RCAL1_height(1, i) < CAL1_position_threshold && RDM2_height(1, i) < DM2_position_threshold && LDM2_height(1, i) > DM2_position_threshold 
%RCAL1_height(1, i) > CAL1_position_threshold && LCAL1_height(1, i) < CAL1_position_threshold && LDM2_height(1, i) < DM2_position_threshold && RDM2_height(1, i) > DM2_position_threshold


for i = 1:time
    if  LCAL1_data_velocity_result(i,1) > velocity_threshold && RCAL1_data_velocity_result(i,1) < velocity_threshold && RDM2_data_velocity_result(i,1) < velocity_threshold && LDM2_data_velocity_result(i,1) > velocity_threshold && LCAL1_height(1, i) > CAL1_position_threshold && RCAL1_height(1, i) < CAL1_position_threshold && RDM2_height(1, i) < DM2_position_threshold && LDM2_height(1, i) > DM2_position_threshold
        FF_BHO(i) = 1;
    elseif  RCAL1_data_velocity_result(i,1) > velocity_threshold && LCAL1_data_velocity_result(i,1) < velocity_threshold && LDM2_data_velocity_result(i,1) < velocity_threshold && RDM2_data_velocity_result(i,1) > velocity_threshold && RCAL1_height(1, i) > CAL1_position_threshold && LCAL1_height(1, i) < CAL1_position_threshold && LDM2_height(1, i) < DM2_position_threshold && RDM2_height(1, i) > DM2_position_threshold
        FF_BHO(i) = 1;
    else
        FF_BHO(i) = 0;
    end
end

%%% Heel off to before toe off  (HO_BTO) %%%
%%% velocity of L(R)CAL1 above 150 mm/s  LCAL1_data_velocity_result(i,1) > 150
%%% velocity of R(L)CAL1 above 150 mm/s  RCAL1_data_velocity_result(i,1) > 150
%%% velocity of R(L)DM2 below 150 mm/s   RDM2_data_velocity_result(i,1) < 150 
%%% velocity of L(R)DM2 above 150 mm/s   LDM2_data_velocity_result(i,1) > 150 
%%% L(R)CAL1 above 35 mm                 LCAL1_height(1, i) > 35 
%%% R(L)CAL1 above 35 mm                 RCAL1_height(1, i) > 35
%%% R(L)DM2 below 35 mm                  RDM2_height(1, i) < 35
%%% L(R)DM2 above 35 mm                  LDM2_height(1, i) > 35

%Velocity rules
%LCAL1_data_velocity_result(i,1) > velocity_threshold && RCAL1_data_velocity_result(i,1) > velocity_threshold && RDM2_data_velocity_result(i,1) < velocity_threshold && LDM2_data_velocity_result(i,1) > velocity_threshold
%RCAL1_data_velocity_result(i,1) > velocity_threshold && LCAL1_data_velocity_result(i,1) > velocity_threshold && LDM2_data_velocity_result(i,1) < velocity_threshold && RDM2_data_velocity_result(i,1) > velocity_threshold

%Position rules
%LCAL1_height(1, i) > CAL1_position_threshold && RCAL1_height(1, i) > CAL1_position_threshold && RDM2_height(1, i) < DM2_position_threshold && LDM2_height(1, i) > DM2_position_threshold 
%RCAL1_height(1, i) > CAL1_position_threshold && LCAL1_height(1, i) > CAL1_position_threshold && LDM2_height(1, i) < DM2_position_threshold && RDM2_height(1, i) > DM2_position_threshold 

for i = 1:time
    if LCAL1_height(1, i) > CAL1_position_threshold && RCAL1_height(1, i) > CAL1_position_threshold && RDM2_height(1, i) < DM2_position_threshold && LDM2_height(1, i) > DM2_position_threshold
        HO_BTO(i) = 1;
    elseif RCAL1_height(1, i) > CAL1_position_threshold && LCAL1_height(1, i) > CAL1_position_threshold && LDM2_height(1, i) < DM2_position_threshold && RDM2_height(1, i) > DM2_position_threshold
        HO_BTO(i) = 1;
    else
        HO_BTO(i) = 0;
    end
end

plot(HS_BFF)
hold on
plot(FF_BHO)
plot(HO_BTO)
legend('HS_BFF', 'FF_BHO', 'HO_BTO')


for i = 1:time
    if LCAL1_height(1, i) < CAL1_position_threshold || RCAL1_height(1, i) < CAL1_position_threshold
        test(i,1) = 1;
    else
        test(i,1) = 0;
    end
    if RCAL1_height(1, i) > CAL1_position_threshold || LCAL1_height(1, i) > CAL1_position_threshold
        test(i,2) = 1;
    else
        test(i,2) = 0;
    end
    if RDM2_height(1, i) < DM2_position_threshold || LDM2_height(1, i) < DM2_position_threshold
        test(i,3) = 1;
    else
        test(i,3) = 0;
    end
    if LDM2_height(1, i) > DM2_position_threshold || RDM2_height(1, i) > DM2_position_threshold
        test(i,4) = 1;
    else
        test(i,4) = 0;
    end
    if LCAL1_data_velocity_result(i,1) > velocity_threshold || RCAL1_data_velocity_result(i,1) > velocity_threshold
        test(i,5) = 1;
    else
        test(i,5) = 0;
    end
    if RCAL1_data_velocity_result(i,1) < velocity_threshold || LCAL1_data_velocity_result(i,1) < velocity_threshold
        test(i,6) = 1;
    else
        test(i,6) = 0;
    end
    if RDM2_data_velocity_result(i,1) < velocity_threshold || LDM2_data_velocity_result(i,1) < velocity_threshold
        test(i,7) = 1;
    else
        test(i,7) = 0;
    end
    if LDM2_data_velocity_result(i,1) > velocity_threshold || RDM2_data_velocity_result(i,1) > velocity_threshold
        test(i,8) = 1;
    else
        test(i,8) = 0;
    end
    if LCAL1_height(1, i) > CAL1_position_threshold || RCAL1_height(1, i) > CAL1_position_threshold
        test(i,9) = 1;
    else
        test(i,9) = 0;
    end
    if RCAL1_height(1, i) < CAL1_position_threshold || LCAL1_height(1, i) < CAL1_position_threshold
        test(i,10) = 1;
    else
        test(i,10) = 0;
    end
    if RDM2_height(1, i) < DM2_position_threshold || LDM2_height(1, i) < DM2_position_threshold
        test(i,11) = 1;
    else
        test(i,11) = 0;
    end
    if LDM2_height(1, i) > DM2_position_threshold || RDM2_height(1, i) > DM2_position_threshold
        test(i,12) = 1;
    else
        test(i,12) = 0;
    end
    if LCAL1_height(1, i) > CAL1_position_threshold || RCAL1_height(1, i) > CAL1_position_threshold
        test(i,13) = 1;
    else
        test(i,13) = 0;
    end
    if RCAL1_height(1, i) > CAL1_position_threshold || LCAL1_height(1, i) > CAL1_position_threshold
        test(i,14) = 1;
    else
        test(i,14) = 0;
    end
    if RDM2_height(1, i) < DM2_position_threshold || LDM2_height(1, i) < DM2_position_threshold
        test(i,15) = 1;
    else
        test(i,15) = 0;
    end
    if LDM2_height(1, i) > DM2_position_threshold || RDM2_height(1, i) > DM2_position_threshold
        test(i,16) = 1;
    else
        test(i,16) = 0;
    end
end

saved_path_name = "C:\Users\a1003\OneDrive\桌面\Thesis\swing_stance_test.csv";
writematrix(test, saved_path_name)