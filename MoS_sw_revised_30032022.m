%%
%%% refer article: Ohtsu 2019 Investigation of balance strategy over gait cycle based on margin of stability
%%% refer article: https://royalsocietypublishing.org/doi/pdf/10.1098/rsif.2020.0194
%%% velocity -> using the window to avoid the jumping 
%%% position unit
%%% check through the video
%% load the data
clear
addpath 'C:\Users\a1003\OneDrive\桌面\Project_Review\side_walk_new\SW\raw_data'
addpath 'C:\Users\a1003\OneDrive\桌面\Project_Review\side_walk_new\SW'
filenames = dir('C:\Users\a1003\OneDrive\桌面\Project_Review\side_walk_new\SW\raw_data');

for n = 1:length(filenames)
    
    try
        subjfile = filenames(n).name;
        subject = subjfile(1:end-4);
        filename = load(filenames(n).name);
        name = filename.(subsref(fieldnames(filename),substruct('{}',{1})));
        label = filename.(subsref(fieldnames(filename),substruct('{}',{1}))).Trajectories.Labeled.Labels;
        path = filename.(subsref(fieldnames(filename),substruct('{}',{1}))).Trajectories.Labeled.Labels;
        sub = subjfile(1:5) + "\";
        % corresponding time
        time = length(name.Trajectories.Labeled.Data(26,1,:));;
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
        CAL1_position_threshold = 35;
        DM2_position_threshold = 40;

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

        end


        C(1,1:4) = [0.9654, 0.9681, 0.9997, 1.0000];
        C(2,1:4) = [ 1.0000, 0, 1.0000, 0];
        C(3,1:4) = [0.6495, 1, 1, 0];

        [~,idx] = pdist2(C,test,'euclidean','Smallest',1);

        t_dummy = zeros([length(idx), 3]);

        for i = 1:length(idx)
            if idx(i) == 1
                t_dummy(i,1) = 1;
            elseif idx(i) == 2
                t_dummy(i,2) = 1;
            elseif idx(i) == 3
                t_dummy(i,3) = 1;
            end
        end
        %{
        figure
        plot(t_dummy(:,1), 'Color', 'blue')
        hold on
        plot(t_dummy(:,2), 'Color', 'red')
        plot(t_dummy(:,3), 'Color', 'green')
        %}

        %legend("Single stance phase (one foot flat)",  "Double stance phase (two foot flat)", "Double stance phase (one foot flat another forefoot on the ground)")
        DLS = t_dummy(:,2);
        SLS = t_dummy(:,1);
        FFO = t_dummy(:,3);


        SS_per = num2str((sum(SLS))/(sum(DLS) + sum(SLS)+ sum(FFO)));
        DS_per = num2str((sum(DLS) + sum(FFO))/(sum(DLS) + sum(SLS)+ sum(FFO)));

        %%% HS_BFF phase, involving markers => R(L)LMAL, R(L)MMAL, L(R)DH, L(R)DM1, L(R)DM5
        RLMAL_position = find(strcmp( path, 'RLMAL'));
        RLMAL_xy_traj = name.Trajectories.Labeled.Data(RLMAL_position,1:2,:);
        RLMAL_xy_traj = reshape(RLMAL_xy_traj, [2, time]);

        RMMAL_position = find(strcmp( path, 'RMMAL'));
        RMMAL_xy_traj = name.Trajectories.Labeled.Data(RMMAL_position,1:2,:);
        RMMAL_xy_traj = reshape(RMMAL_xy_traj, [2, time]);

        LLMAL_position = find(strcmp( path, 'LLMAL'));
        LLMAL_xy_traj = name.Trajectories.Labeled.Data(LLMAL_position,1:2,:);
        LLMAL_xy_traj = reshape(LLMAL_xy_traj, [2, time]);

        LMMAL_position = find(strcmp( path, 'LMMAL'));
        LMMAL_xy_traj = name.Trajectories.Labeled.Data(LMMAL_position,1:2,:);
        LMMAL_xy_traj = reshape(LMMAL_xy_traj, [2, time]);

        RDH_position = find(strcmp( path, 'RDH'));
        RDH_xy_traj = name.Trajectories.Labeled.Data(RDH_position,1:2,:);
        RDH_xy_traj = reshape(RDH_xy_traj, [2, time]);

        LDH_position = find(strcmp( path, 'LDH'));
        LDH_xy_traj = name.Trajectories.Labeled.Data(LDH_position,1:2,:);
        LDH_xy_traj = reshape(LDH_xy_traj, [2, time]);

        RDM1_position = find(strcmp( path, 'RDM1'));
        RDM1_xy_traj = name.Trajectories.Labeled.Data(RDM1_position,1:2,:);
        RDM1_xy_traj = reshape(RDM1_xy_traj, [2, time]);

        LDM1_position = find(strcmp( path, 'LDM1'));
        LDM1_xy_traj = name.Trajectories.Labeled.Data(LDM1_position,1:2,:);
        LDM1_xy_traj = reshape(LDM1_xy_traj, [2, time]);

        RPM1_position = find(strcmp( path, 'RPM1'));
        RPM1_xy_traj = name.Trajectories.Labeled.Data(RPM1_position,1:2,:);
        RPM1_xy_traj = reshape(RPM1_xy_traj, [2, time]);

        LPM1_position = find(strcmp( path, 'LPM1'));
        LPM1_xy_traj = name.Trajectories.Labeled.Data(LPM1_position,1:2,:);
        LPM1_xy_traj = reshape(LPM1_xy_traj, [2, time]);

        RDM5_position = find(strcmp( path, 'RDM5'));
        RDM5_xy_traj = name.Trajectories.Labeled.Data(RDM5_position,1:2,:);
        RDM5_xy_traj = reshape(RDM5_xy_traj, [2, time]);

        LDM5_position = find(strcmp( path, 'LDM5'));
        LDM5_xy_traj = name.Trajectories.Labeled.Data(LDM5_position,1:2,:);
        LDM5_xy_traj = reshape(LDM5_xy_traj, [2, time]);

        RPM5_position = find(strcmp( path, 'RPM5'));
        RPM5_xy_traj = name.Trajectories.Labeled.Data(RPM5_position,1:2,:);
        RPM5_xy_traj = reshape(RPM5_xy_traj, [2, time]);

        LPM5_position = find(strcmp( path, 'LPM5'));
        LPM5_xy_traj = name.Trajectories.Labeled.Data(LPM5_position,1:2,:);
        LPM5_xy_traj = reshape(LPM5_xy_traj, [2, time]);



        %% final step of the calculation of MoS
        %%%% decision of the phase => compare the distance => decide the value

        %%%%%%%%%%%  DLS  %%%%%%%%%%%%%

        for i = 1:time
            %%%%% both foot on the ground
            if DLS(i) == 1
                xcom = XCOM*1000;
                mid_forefoot = (RDM5_xy_traj + LDM5_xy_traj)/2;
                mid_heel = (RLMAL_xy_traj + LLMAL_xy_traj)/2;
                AP_d_1 = abs(xcom(2,i) - mid_forefoot(2,i));
                AP_d_2 = abs(xcom(2,i) - mid_heel(2,i));
                AP_d_3 = abs(mid_forefoot(2,i) - mid_heel(2,i));
                if (AP_d_1 + AP_d_2) == AP_d_3
                    AP_MoS = min([AP_d_1, AP_d_2]);
                    AP_MoS_DLS(i) = AP_MoS;
                else
                    AP_MoS = min([AP_d_1, AP_d_2]);
                    AP_MoS_DLS(i) = 0-AP_MoS;
                end

                mid_right = (RDM5_xy_traj + RLMAL_xy_traj)/2;
                mid_left = (LDM5_xy_traj + LLMAL_xy_traj)/2;
                ML_d_1 = abs(xcom(1,i) - mid_right(1,i));
                ML_d_2 = abs(xcom(1,i) - mid_left(1,i));
                ML_d_3 = abs(mid_right(1,i) - mid_left(1,i));
                if (ML_d_1 + ML_d_2) == ML_d_3
                    ML_MoS = min([ML_d_1, ML_d_2]);
                    ML_MoS_DLS(i) = ML_MoS;
                else
                    ML_MoS = min([ML_d_1, ML_d_2]);
                    ML_MoS_DLS(i) = 0-ML_MoS;
                end
            else 
                AP_MoS_DLS(i) = 0;
                ML_MoS_DLS(i) = 0;
            end
        end


        %%%%% both forefoot on the ground and one heel in the air
        for i = 1:time
            %%%% right heel off
            if  FFO(i) == 1 && RCAL1_height(i)> LCAL1_height(i)
                xcom = XCOM*1000;
                mid_forefoot = (RDM5_xy_traj + LDM5_xy_traj)/2;
                mid_heel = (RPM5_xy_traj + LLMAL_xy_traj)/2;
                AP_d_1 = abs(xcom(2,i) - mid_forefoot(2,i));
                AP_d_2 = abs(xcom(2,i) - mid_heel(2,i));
                AP_d_3 = abs(mid_forefoot(2,i) - mid_heel(2,i));
                if (AP_d_1 + AP_d_2) == AP_d_3
                    AP_MoS = min([AP_d_1, AP_d_2]);
                    AP_MoS_FFO(i) = AP_MoS;
                else
                    AP_MoS = min([AP_d_1, AP_d_2]);
                    AP_MoS_FFO(i) = 0-AP_MoS;
                end

                mid_right = (RDM5_xy_traj + RPM5_xy_traj)/2;
                mid_left = (LDM5_xy_traj + LLMAL_xy_traj)/2;
                ML_d_1 = abs(xcom(1,i) - mid_right(1,i));
                ML_d_2 = abs(xcom(1,i) - mid_left(1,i));
                ML_d_3 = abs(mid_right(1,i) - mid_left(1,i));
                if (ML_d_1 + ML_d_2) == ML_d_3
                    ML_MoS = min([ML_d_1, ML_d_2]);
                    ML_MoS_FFO(i) = ML_MoS;
                else
                    ML_MoS = min([ML_d_1, ML_d_2]);
                    ML_MoS_FFO(i) = 0-ML_MoS;
                end
            %%%% left heel off
            elseif  FFO(i) == 1 && LCAL1_height(i)> RCAL1_height(i)
                xcom = XCOM*1000;
                mid_forefoot = (LDM5_xy_traj + RDM5_xy_traj)/2;
                mid_heel = (LPM5_xy_traj + RLMAL_xy_traj)/2;
                AP_d_1 = abs(xcom(2,i) - mid_forefoot(2,i));
                AP_d_2 = abs(xcom(2,i) - mid_heel(2,i));
                AP_d_3 = abs(mid_forefoot(2,i) - mid_heel(2,i));
                if (AP_d_1 + AP_d_2) == AP_d_3
                    AP_MoS = min([AP_d_1, AP_d_2]);
                    AP_MoS_FFO(i) = AP_MoS;
                else
                    AP_MoS = min([AP_d_1, AP_d_2]);
                    AP_MoS_FFO(i) = 0-AP_MoS;
                end

                mid_right = (RDM5_xy_traj + RLMAL_xy_traj )/2;
                mid_left = (LDM5_xy_traj + LPM5_xy_traj)/2;
                ML_d_1 = abs(xcom(1,i) - mid_right(1,i));
                ML_d_2 = abs(xcom(1,i) - mid_left(1,i));
                ML_d_3 = abs(mid_right(1,i) - mid_left(1,i));
                if (ML_d_1 + ML_d_2) == ML_d_3
                    ML_MoS = min([ML_d_1, ML_d_2]);
                    ML_MoS_FFO(i) = ML_MoS;
                else
                    ML_MoS = min([ML_d_1, ML_d_2]);
                    ML_MoS_FFO(i) = 0-ML_MoS;
                end    
            else
                AP_MoS_FFO(i) = 0;
                ML_MoS_FFO(i) = 0;
            end
        end
        %%%%% single limb support
        for i = 1:time
            %%%% right foot off
            if  SLS(i) == 1 && RCAL1_height(i)> LCAL1_height(i)
                xcom = XCOM*1000;
                mid_forefoot = (LDM1_xy_traj + LDM5_xy_traj)/2;
                mid_heel = (LMMAL_xy_traj + LLMAL_xy_traj)/2;
                AP_d_1 = abs(xcom(2,i) - mid_forefoot(2,i));
                AP_d_2 = abs(xcom(2,i) - mid_heel(2,i));
                AP_d_3 = abs(mid_forefoot(2,i) - mid_heel(2,i));
                if (AP_d_1 + AP_d_2) == AP_d_3
                    AP_MoS = min([AP_d_1, AP_d_2]);
                    AP_MoS_SLS(i) = AP_MoS;
                else
                    AP_MoS = min([AP_d_1, AP_d_2]);
                    AP_MoS_SLS(i) = 0-AP_MoS;
                end

                mid_right = (LDM1_xy_traj + LMMAL_xy_traj)/2;
                mid_left = (LDM5_xy_traj + LLMAL_xy_traj)/2;
                ML_d_1 = abs(xcom(1,i) - mid_right(1,i));
                ML_d_2 = abs(xcom(1,i) - mid_left(1,i));
                ML_d_3 = abs(mid_right(1,i) - mid_left(1,i));
                if (ML_d_1 + ML_d_2) == ML_d_3
                    ML_MoS = min([ML_d_1, ML_d_2]);
                    ML_MoS_SLS(i) = ML_MoS;
                else
                    ML_MoS = min([ML_d_1, ML_d_2]);
                    ML_MoS_SLS(i) = 0-ML_MoS;
                end
            %%%% left foot off
            elseif  SLS(i) == 1 && LCAL1_height(i)> RCAL1_height(i)
                xcom = XCOM*1000;
                mid_forefoot = (RDM1_xy_traj + RDM5_xy_traj)/2;
                mid_heel = (RMMAL_xy_traj + RLMAL_xy_traj)/2;
                AP_d_1 = abs(xcom(2,i) - mid_forefoot(2,i));
                AP_d_2 = abs(xcom(2,i) - mid_heel(2,i));
                AP_d_3 = abs(mid_forefoot(2,i) - mid_heel(2,i));
                if (AP_d_1 + AP_d_2) == AP_d_3
                    AP_MoS = min([AP_d_1, AP_d_2]);
                    AP_MoS_SLS(i) = AP_MoS;
                else
                    AP_MoS = min([AP_d_1, AP_d_2]);
                    AP_MoS_SLS(i) = 0-AP_MoS;
                end

                mid_right = (RDM5_xy_traj + RLMAL_xy_traj )/2;
                mid_left = (RDM1_xy_traj + RMMAL_xy_traj)/2;
                ML_d_1 = abs(xcom(1,i) - mid_right(1,i));
                ML_d_2 = abs(xcom(1,i) - mid_left(1,i));
                ML_d_3 = abs(mid_right(1,i) - mid_left(1,i));
                if (ML_d_1 + ML_d_2) == ML_d_3
                    ML_MoS = min([ML_d_1, ML_d_2]);
                    ML_MoS_SLS(i) = ML_MoS;
                else
                    ML_MoS = min([ML_d_1, ML_d_2]);
                    ML_MoS_SLS(i) = 0-ML_MoS;
                end    
            else
                AP_MoS_SLS(i) = 0;
                ML_MoS_SLS(i) = 0;
            end
        end


        AP_MoS_DS = NaN(1,time);
        AP_MoS_SS = NaN(1,time);
        for i = 1:time
            if AP_MoS_DLS(1,i) ~= 0
                AP_MoS_all(i) = AP_MoS_DLS(1,i);
                AP_MoS_DS(1, i) = AP_MoS_DLS(1,i);
            elseif AP_MoS_SLS(1,i) ~= 0
                AP_MoS_all(i) = AP_MoS_SLS(1,i);
                AP_MoS_SS(1, i) = AP_MoS_SLS(1,i);
            elseif AP_MoS_FFO(1,i) ~= 0
                AP_MoS_all(i) = AP_MoS_FFO(1,i);
                AP_MoS_DS(1, i) = AP_MoS_FFO(1,i);
            end
        end

        [val, pos] = findpeaks(0-AP_MoS_all(1:end-10), "MinPeakProminence", 150, "MinPeakDistance", 50);
        Min_Mean_AP = num2str(0-mean(val));
        DS_MoS_Mean_AP = num2str(nanmean(AP_MoS_DS(1:end-10)));
        DS_MoS_Var_AP = num2str(nanstd(AP_MoS_DS(1:end-10)));
        SS_MoS_Mean_AP = num2str(nanmean(AP_MoS_SS(1:end-10)));
        SS_MoS_Var_AP = num2str(nanstd(AP_MoS_SS(1:end-10)));

        ML_MoS_DS = NaN(1,time);
        ML_MoS_SS = NaN(1,time);
        for i = 1:time
            if ML_MoS_DLS(1,i) ~= 0
                ML_MoS_all(i) = ML_MoS_DLS(1,i);
                ML_MoS_DS(1, i) = ML_MoS_DLS(1,i);
            elseif ML_MoS_SLS(1,i) ~= 0
                ML_MoS_all(i) = ML_MoS_SLS(1,i);
                ML_MoS_SS(1, i) = ML_MoS_SLS(1,i);
            elseif ML_MoS_FFO(1,i) ~= 0
                ML_MoS_all(i) = ML_MoS_FFO(1,i);
                ML_MoS_DS(1, i) = ML_MoS_FFO(1,i);
            end
        end
        [val, pos] = findpeaks(0-ML_MoS_all(1:end-10), "MinPeakProminence", 150, "MinPeakDistance", 50);
        Min_Mean_ML = num2str(0-mean(val));
        DS_MoS_Mean_ML = num2str(nanmean(ML_MoS_DS(1:end-10)));
        DS_MoS_Var_ML = num2str(nanstd(ML_MoS_DS(1:end-10)));
        SS_MoS_Mean_ML = num2str(nanmean(ML_MoS_SS(1:end-10)));
        SS_MoS_Var_ML = num2str(nanstd(ML_MoS_SS(1:end-10)));


        col = ["subject", "Single stance", "Double stance", "Min_Mean_AP", "DS_MoS_Mean_AP", ...
            "DS_MoS_Var_AP", "SS_MoS_Mean_AP", "SS_MoS_Var_AP", "Min_Mean_ML", "DS_MoS_Mean_ML", ...
            "DS_MoS_Var_ML", "SS_MoS_Mean_ML", "SS_MoS_Var_ML"];
        infor = [string(subject(1:3)), string(SS_per), string(DS_per), string(Min_Mean_AP), string(DS_MoS_Mean_AP), ...
            string(DS_MoS_Var_AP), string(SS_MoS_Mean_AP), string(SS_MoS_Var_AP), string(Min_Mean_ML), string(DS_MoS_Mean_ML), ...
            string(DS_MoS_Var_ML), string(SS_MoS_Mean_ML), string(SS_MoS_Var_ML)];
        save_data = [col;infor];

        saved_path_name = "C:/Users/a1003/OneDrive/桌面/Project_Review/side_walk_new/SW/csv_file/MoS/SW/" + subject + "_MoS_SW.csv";
        writematrix(save_data, saved_path_name)
        clear
        addpath 'C:\Users\a1003\OneDrive\桌面\Project_Review\side_walk_new\SW\raw_data'
        addpath 'C:\Users\a1003\OneDrive\桌面\Project_Review\side_walk_new\SW'
        filenames = dir('C:\Users\a1003\OneDrive\桌面\Project_Review\side_walk_new\SW\raw_data');   


        %{
        figure
        plot(AP_MoS_DLS/10)
        hold on
        plot(AP_MoS_FFO/10)
        plot(AP_MoS_SLS/10)
        legend('AP_MoS_DLS', 'AP_MoS_FFO', 'AP_MoS_SLS')
        xlabel("time")
        ylabel("MoS (cm)")


        figure
        plot(ML_MoS_DLS/10)
        hold on
        plot(ML_MoS_FFO/10)
        plot(ML_MoS_SLS/10)
        legend('ML_MoS_DLS', 'ML_MoS_FFO', 'ML_MoS_SLS')
        xlabel("time")
        ylabel("MoS (cm)")
        %}
    catch
        subjfile
        %clear
        addpath 'C:\Users\a1003\OneDrive\桌面\Project_Review\side_walk_new\SW\raw_data'
        addpath 'C:\Users\a1003\OneDrive\桌面\Project_Review\side_walk_new\SW'
        filenames = dir('C:\Users\a1003\OneDrive\桌面\Project_Review\side_walk_new\SW\raw_data'); 
        continue
    end
end        
