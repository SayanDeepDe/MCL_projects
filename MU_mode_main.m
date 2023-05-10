%% This code is used to compute MU-synergies (intra muscle) stabilizing total force during drifts
clc
clear
close all
subjid='subj16';
p=cd(subjid);

%% Load data
load drift_data.mat
L_drift=length(data);
drift_f=cell(length(data),1);
for i=1:length(data)
    drift_f{i}.signals=data{i}.signals;
    drift_f{i}.name=data{i}.name;
end

load MVC_data.mat
L_MVC=length(data);
MVC=cell(length(data),1);
for i=1:L_MVC
    MVC{i}.signals=data{i}.signals;
    MVC{i}.name=data{i}.name;
end

load enslaving_data.mat
L_enslaving=length(data);
ensl=cell(length(data),1);
for i=1:L_enslaving
    ensl{i}.signals=data{i}.signals;
    ensl{i}.name=data{i}.name;
end

load sensor1.mat
L_sensor1=length(sensor1_data);
sensor1=cell(L_sensor1,1);
for i=1:L_sensor1
    sensor1{i}.signals=sensor1_data{i}.signals;
    sensor1{i}.name=sensor1_data{i}.name;
end

load sensor2.mat
L_sensor2=length(sensor2_data);
sensor2=cell(L_sensor2,1);
for i=1:L_sensor2
    sensor2{i}.signals=sensor2_data{i}.signals;
    sensor2{i}.name=sensor2_data{i}.name;
end

%% initializing parameters
critical_n_MU=5; % critical number of MU detected to accept the trials
n_drift=5;       % number of episodes of drift which is pre-fixed before data collection
%% matching the Force, sensor1 and sensor2 data. Contain only the trials which has both sensor1 and sensor2
[match_F_sensor,match_sensor1,match_sensor2,null_sensor1,null_sensor2,L_match]= match(drift_f,sensor1,sensor2,L_drift,L_sensor1,L_sensor2,critical_n_MU);

%% Enslaving matrix
ENSL_L=nan(4,4); %initializing the ENSL matrix with nan for lt and rt hand
ENSL_R=nan(4,4);
 for i=1:length(ensl)
     [k_I,k_M,k_R,k_L,lt,rt,n]=enslave_coeff(ensl{i}.signals,data{i}.name);
     if lt==1 %building the enslaving matrix based on task finger (n=1(I), 2(M), 3(R), 4(L))
         ENSL_L(1,n)=k_I;
         ENSL_L(2,n)=k_M;
         ENSL_L(3,n)=k_R;
         ENSL_L(4,n)=k_L;
     elseif rt==1
         ENSL_R(1,n)=k_I;
         ENSL_R(2,n)=k_M;
         ENSL_R(3,n)=k_R;
         ENSL_R(4,n)=k_L;

     end

 end

%% computing the mode data. m = E^(-1) * F
 M_data=cell(L_enslaving,1);
 for i=1:L_enslaving
     M_data{i}.name=ensl{i}.name;
     if contains(M_data{i}.name,"Lt","IgnoreCase",true)
        M_data{i}.signals=(inv(ENSL_L) * (ensl{i}.signals(:,2:5))')';
     else
        M_data{i}.signals=(inv(ENSL_R) * (ensl{i}.signals(:,2:5))')'; 
     end
     
 end

%% split = splitting the episodes
t_drift=cell(L_drift,1);
f_drift=cell(L_drift,1);
for i=1:L_drift
    t_drift{i}.signals=match_F_sensor{i}.signals(:,1);
    t_drift{i}.name=match_F_sensor{i}.name;
    f_drift{i}.signals=match_F_sensor{i}.signals(:,2:5);
    f_drift{i}.name=match_F_sensor{i}.name;
end

%% Entering each trial and processing them individually

sensor1_cycle=cell(L_match,1);          % store cyclical data for each trial from sensor1
sensor2_cycle=cell(L_match,1);          % store cyclical data for each trial from sensor2
T_cycle=cell(L_match,1);                % store cyclical time points for each trial from raw Force data
F_cycle=cell(L_match,1);                % store cyclical force data for each trial from raw Force data
F_cycle_tot=cell(L_match,1);            % store cyclical sum of 4-finger forces for each trial from Force data
T_drift=cell(L_match,n_drift);          % store time points of drift episodes for each trial from Force data
F_drift=cell(L_match,n_drift);          % store force of drift episodes for each trial from Force data
FT_drift=cell(L_match,n_drift);
T_sensor_drift=cell(L_match,n_drift);   % store time points of drift episodes for each trial from EMG file
sensor1_drift=cell(L_match,n_drift);    % store f-MU of drift episodes for each trial from sensor1
mode_sensor1_drift=cell(L_match,n_drift); % store Mode data of drift episodes for each trial from sensor1
sensor1_T_drift=cell(L_match,n_drift);
sensor2_drift=cell(L_match,n_drift);    % store f-MU of drift episodes for each trial from sensor2
mode_sensor2_drift=cell(L_match,n_drift); % store Mode data of drift episodes for each trial from sensor2
sensor2_T_drift=cell(L_match,n_drift);
mode_sensor12_drift=cell(L_match,n_drift); % store Mode data of drift episodes for each trial from sensor1+2



for j=1:L_match % each trial
   [sensor1_cyc,sensor2_cyc,T_sensor_cyc,sensor1_cyc_ave,sensor1_cyc_ave_new,sensor2_cyc_ave,sensor2_cyc_ave_new,T_cyc_ave,T_sensor_cyc_ave,F_cyc_ave,F_cyc_tot_ave,F_cyc_tot_ave_new,T_sensor_d,sensor1_d,sensor2_d,F_d,T_d,FT_d,sensor1_T_d,sensor2_T_d]=split_drift(match_sensor1{j}.signals,match_sensor2{j}.signals,match_F_sensor{j}.signals,match_F_sensor{j}.name,ENSL_L,ENSL_R,MVC);

%   For PCA and jacobian, we use only the cyclical force and EMG.
%   PCA and Jacobian for sensor1
    [coeff_sensor1,score,latent,tsquared,explained] = pca(sensor1_cyc_ave); %pca
    mode_sensor1=(sensor1_cyc-mean(sensor1_cyc))*coeff_sensor1(:,1:2);  %mode magnitude
    mode_sensor1_dt=detrend(mode_sensor1); %detrend the mode magnitude
    mode_sensor1_ave=(sensor1_cyc_ave-mean(sensor1_cyc_ave))*coeff_sensor1(:,1:2);
    jacob_sensor1(j,:)=regress(diff(F_cyc_tot_ave),diff(mode_sensor1_ave));       %jacobian for mode magnitude    
    jacob_sensor1_VAF_old(j,:)=fitlm(diff(mode_sensor1_ave), diff(F_cyc_tot_ave)).Rsquared.Adjusted;
%   Jacobian for sensor1 NEW
    mode_sensor1_ave_new=(sensor1_cyc_ave_new-mean(sensor1_cyc_ave_new))*coeff_sensor1(:,1:2);
    jacob_sensor1_new(j,:)=regress(diff(F_cyc_tot_ave_new),diff(mode_sensor1_ave_new));       %jacobian for mode magnitude    
    jacob_sensor1_VAF_new(j,:)=fitlm(diff(mode_sensor1_ave_new), diff(F_cyc_tot_ave_new)).Rsquared.Adjusted;

%   PCA and Jacobian for sensor2
    [coeff_sensor2,score,latent,tsquared,explained] = pca(sensor2_cyc_ave); %pca
    mode_sensor2=(sensor2_cyc-mean(sensor2_cyc))*coeff_sensor2(:,1:2);
    mode_sensor2_dt=detrend(mode_sensor2);                              %detrend the mode magnitude
    mode_sensor2_ave=(sensor2_cyc_ave-mean(sensor2_cyc_ave))*coeff_sensor2(:,1:2);  %mode magnitude
    jacob_sensor2(j,:)=regress(diff(F_cyc_tot_ave),diff(mode_sensor2_ave));       %jacobian for mode magnitude
    jacob_sensor2_VAF_old(j,:)=fitlm(diff(mode_sensor2_ave), diff(F_cyc_tot_ave)).Rsquared.Adjusted;
%   Jacobian for sensor2 NEW
    mode_sensor2_ave_new=(sensor2_cyc_ave_new-mean(sensor2_cyc_ave_new))*coeff_sensor2(:,1:2);
    jacob_sensor2_new(j,:)=regress(diff(F_cyc_tot_ave_new),diff(mode_sensor2_ave_new));       %jacobian for mode magnitude    
    jacob_sensor2_VAF_new(j,:)=fitlm(diff(mode_sensor2_ave_new), diff(F_cyc_tot_ave_new)).Rsquared.Adjusted;

%   PCA and Jacobian for sensor1+2
    sensor12_cyc_ave=cat(2,sensor1_cyc_ave,sensor2_cyc_ave);
    sensor12_cyc=cat(2,sensor1_cyc,sensor2_cyc);
    [coeff_sensor12,score,latent,tsquared,explained] = pca(sensor12_cyc_ave); %pca
    mode_sensor12=(sensor12_cyc-mean(sensor12_cyc))*coeff_sensor12(:,1:2);  %mode magnitude
    mode_sensor12_dt=detrend(mode_sensor12);                              %detrend the mode magnitude
    mode_sensor12_ave=(sensor12_cyc_ave-mean(sensor12_cyc_ave))*coeff_sensor12(:,1:2);
    jacob_sensor12(j,:)=regress(diff(F_cyc_tot_ave),diff(mode_sensor12_ave));       %jacobian for mode magnitude
    jacob_sensor12_VAF_old(j,:)=fitlm(diff(mode_sensor12_ave), diff(F_cyc_tot_ave)).Rsquared.Adjusted;
%   Jacobian for sensor1+2 NEW
    sensor12_cyc_ave_new=cat(2,sensor1_cyc_ave_new,sensor2_cyc_ave_new);
    mode_sensor12_ave_new=(sensor12_cyc_ave_new-mean(sensor12_cyc_ave_new))*coeff_sensor12(:,1:2);
    jacob_sensor12_new(j,:)=regress(diff(F_cyc_tot_ave_new),diff(mode_sensor12_ave_new));       %jacobian for mode magnitude    
    jacob_sensor12_VAF_new(j,:)=fitlm(diff(mode_sensor12_ave_new), diff(F_cyc_tot_ave_new)).Rsquared.Adjusted;

%     figure        % plot for loading factors
%     bar(coeff_sensor12(:,1:2),'DisplayName','coeff_sensor12(:,1:2)');
%     xline(size(sensor1_cyc_ave,2));

%   storing chopped drift trials. Each row=trial, column=episode
    [F_drift{j,:}]=F_d{:};                  %storing Force episodes from each trial
    [T_drift{j,:}]=T_d{:};                  %storing Time(Force) episodes from each trial
    [FT_drift{j,:}]=FT_d{:};
    [T_sensor_drift{j,:}]=T_sensor_d{:};    %storing Time(EMG) episodes from each trial
    [sensor1_drift{j,:}]=sensor1_d{:};      %storing sensor1(EMG) episodes from each trial
    [sensor1_T_drift{j,:}]=sensor1_T_d{:};
    [sensor2_drift{j,:}]=sensor2_d{:};      %storing sensor2(EMG) episodes from each trial
    [sensor2_T_drift{j,:}]=sensor2_T_d{:};

    for i=1:size(sensor1_d,2) %mode data for EMG during drift episodes
    [mode_sensor1_drift{j,i}]=(sensor1_d{:,i})*coeff_sensor1(:,1:2);
    [mode_sensor2_drift{j,i}]=(sensor2_d{:,i})*coeff_sensor2(:,1:2);
    [mode_sensor12_drift{j,i}]=(cat(2,sensor1_d{:,i},sensor2_d{:,i}))*coeff_sensor12(:,1:2);
    end

% figure
% plot(cat(1,T_sensor_drift{j,:}),cat(1,mode_sensor1_drift{j,:}),"r--",cat(1,T_sensor_drift{j,:}),cat(1,mode_sensor2_drift{j,:}),"b:",cat(1,T_sensor_drift{j,:}),cat(1,mode_sensor12_drift{j,:}),"k");

%     figure
%     subplot(3,1,1);
%     plot(T_sensor_cyc,mode_sensor1);
%     title('Mode magnitude')
%     subtitle('Sensor1')
%     subplot(3,1,2);
%     plot(T_sensor_cyc,mode_sensor2);
%     title('Mode magnitude')
%     subtitle('Sensor2')
%     subplot(3,1,3);
%     plot(T_sensor_cyc,mode_sensor12);
%     title('Mode magnitude')
%     subtitle('Sensor1+2')

    sensor1_cycle{j}.signals=sensor1_cyc_ave;
    sensor2_cycle{j}.signals=sensor2_cyc_ave;
    T_cycle{j}.signals=T_cyc_ave;

end


%% converting Force data in drift episodes to Force-mode


 for i=1:size(F_drift,1)
     
     if contains(match_F_sensor{i}.name,"Lt","IgnoreCase",true)
         for j=1:size(F_drift,2)
                 [F_datamode{i,j}]=deal((inv(ENSL_L)* (F_drift{i,j})')');
         end
     else
         for j=1:size(F_drift,2)
                [F_datamode{i,j}]=deal((inv(ENSL_R)* (F_drift{i,j})')'); 
         end
     end
     
 end


%% chopping off 250ms windows pre and post drift



a=[1,2,3];
b=[1,2,3];
c=combvec(a,b);
pre_F_mode=nan(n_drift,4,size(F_drift,1));
post_F_mode=nan(n_drift,4,size(F_drift,1));

for k=1:size(F_drift,1) % moving through each trial.
    pre_F=26750;    % start of drift episode 1
    post_F=43000;   % end of drift episode 1

    for l=1:size(F_drift,2) % moving through each episode of drift of a single trial.
        % start indexes of 3 time windows pre drift
        pre_F_1=pre_F;
        post_F_1=post_F-1250;
       for w=1:3 % moving through each time window
           [~,pre_T_F(1,w)]=min(abs(T_drift{k,l}-pre_F_1)); 
           [~,pre_T_F(2,w)]=min(abs(T_drift{k,l}-(pre_F_1+250)));

           [~,post_T_F(1,w)]=min(abs(T_drift{k,l}-(post_F_1)));
           [~,post_T_F(2,w)]=min(abs(T_drift{k,l}-(post_F_1+250)));

           [~,pre_T_E(1,w)]=min(abs(T_sensor_drift{k,l}-(pre_F_1/1000)));
           [~,pre_T_E(2,w)]=min(abs(T_sensor_drift{k,l}-((pre_F_1+250)/1000)));

           [~,post_T_E(1,w)]=min(abs(T_sensor_drift{k,l}-(post_F_1/1000)));
           [~,post_T_E(2,w)]=min(abs(T_sensor_drift{k,l}-((post_F_1+250)/1000)));

           
           %sensor 1
           pre_E1(w,:)=mean((mode_sensor1_drift{k,l}(pre_T_E(1,w):pre_T_E(2,w),:)));
           post_E1(w,:)=mean((mode_sensor1_drift{k,l}(post_T_E(1,w):post_T_E(2,w),:)));
           pre_M1_sensor1(w,l,k)=pre_E1(w,1);
           pre_M2_sensor1(w,l,k)=pre_E1(w,2);
           post_M1_sensor1(w,l,k)=post_E1(w,1);
           post_M2_sensor1(w,l,k)=post_E1(w,2);


           %sensor 2
           pre_E2(w,:)=mean((mode_sensor2_drift{k,l}(pre_T_E(1,w):pre_T_E(2,w),:)));
           post_E2(w,:)=mean((mode_sensor2_drift{k,l}(post_T_E(1,w):post_T_E(2,w),:)));
           pre_M1_sensor2(w,l,k)=pre_E2(w,1);
           pre_M2_sensor2(w,l,k)=pre_E2(w,2);
           post_M1_sensor2(w,l,k)=post_E2(w,1);
           post_M2_sensor2(w,l,k)=post_E2(w,2);

           %sensor 12
           pre_E12(w,:)=mean((mode_sensor12_drift{k,l}(pre_T_E(1,w):pre_T_E(2,w),:)));
           post_E12(w,:)=mean((mode_sensor12_drift{k,l}(post_T_E(1,w):post_T_E(2,w),:)));
           pre_M1_sensor12(w,l,k)=pre_E12(w,1);
           pre_M2_sensor12(w,l,k)=pre_E12(w,2);
           post_M1_sensor12(w,l,k)=post_E12(w,1);
           post_M2_sensor12(w,l,k)=post_E12(w,2);

           pre_F_1=pre_F_1+500;
           post_F_1=post_F_1+500;
       end

       %Forcemode
           pre_F_mode(l,:,k)=mean((F_datamode{k,l}(pre_T_F(1,3):pre_T_F(2,3),:)));
           post_F_mode(l,:,k)=mean((F_datamode{k,l}(post_T_F(1,3):post_T_F(2,3),:)));

        

        for n=1:size(c,2)   %storing the 9 combinations of difference between 3 post and 3 pre time-windows
            h=1;
            d_M1_sensor1(n,l,k)=post_E1(c(2,n),h)-pre_E1(c(1,n),h);     % rows = 9 combinations of difference between post and pre drift, col = episode number 
            d_M2_sensor1(n,l,k)=post_E1(c(2,n),h+1)-pre_E1(c(1,n),h+1);

            d_M1_sensor2(n,l,k)=post_E2(c(2,n),h)-pre_E2(c(1,n),h);     % rows = 9 combinations of difference between post and pre drift, col = episode number 
            d_M2_sensor2(n,l,k)=post_E2(c(2,n),h+1)-pre_E2(c(1,n),h+1);

            d_M1_sensor12(n,l,k)=post_E12(c(2,n),h)-pre_E12(c(1,n),h);     % rows = 9 combinations of difference between post and pre drift, col = episode number 
            d_M2_sensor12(n,l,k)=post_E12(c(2,n),h+1)-pre_E12(c(1,n),h+1);
        end
        
        F_total_drift{k,l}=sum(F_drift{k,l},2);
        
        pre_F=pre_F+23000;
        post_F=post_F+23000;

    end

     % Horizontally concatenate, pad with NaNs
        maxNumRows = max(cellfun(@(c) length(c), F_total_drift(k,:)));  % max number of columns
        F_total_drift_mat = cell2mat(cellfun(@(c){padarray(c,[maxNumRows-length(c),0],NaN,'Post')}, F_total_drift(k,:)));
        F_datamode_mat = cell2mat(cellfun(@(c){padarray(c,[maxNumRows-length(c),0],NaN,'Post')}, F_datamode(k,:)));
        
        % define the column indices for individual modes (evenly spaced)
        start_index_mode1 = 1; % start indices different for each mode
        start_index_mode2 = 2;
        start_index_mode3 = 3;
        start_index_mode4 = 4;
        end_index = 20;         % same end index for all modes
        step_size = 4;          % same step size for all modes
        group_indices_mode1 = start_index_mode1:step_size:end_index;
        group_indices_mode2 = start_index_mode2:step_size:end_index;
        group_indices_mode3 = start_index_mode3:step_size:end_index;
        group_indices_mode4 = start_index_mode4:step_size:end_index;
        
        % create a new matrix to store the grouped columns
        F_mode1 = nan(size(F_datamode_mat, 1), numel(group_indices_mode1));
        F_mode2 = nan(size(F_datamode_mat, 1), numel(group_indices_mode2));
        F_mode3 = nan(size(F_datamode_mat, 1), numel(group_indices_mode3));
        F_mode4 = nan(size(F_datamode_mat, 1), numel(group_indices_mode4));
        
% loop over the group indices for each mode and extract the corresponding columns from F_datamode.
        for i = 1:numel(group_indices_mode1)
            F_mode1(:, i) = F_datamode_mat(:, group_indices_mode1(i));
            F_mode2(:, i) = F_datamode_mat(:, group_indices_mode2(i));
            F_mode3(:, i) = F_datamode_mat(:, group_indices_mode3(i));
            F_mode4(:, i) = F_datamode_mat(:, group_indices_mode4(i));
        end

        F_total_drift_mean = mean(F_total_drift_mat,2,'omitnan');   % mean of total force across 5 episodes in 1 trial
        F_total_drift_std=std(F_total_drift_mat,0,2,"omitnan");     % std of total force across 5 episodes in 1 trial
        [F_total_drift_pad{k,1}]=F_total_drift_mat;
        [F_total_drift_pad_mean{k,1}]=F_total_drift_mean;
        [F_total_drift_pad_std{k,1}]=F_total_drift_std;

        F_mode1_mean = mean(F_mode1,2,'omitnan');   % mean of force mode 1 across 5 episodes in 1 trial
        F_mode1_std=std(F_mode1,0,2,"omitnan");     % std of force mode 1 across 5 episodes in 1 trial

        F_mode2_mean = mean(F_mode2,2,'omitnan');   % mean of force mode 2 across 5 episodes in 1 trial
        F_mode2_std=std(F_mode2,0,2,"omitnan");     % std of force mode 2 across 5 episodes in 1 trial

        F_mode3_mean = mean(F_mode3,2,'omitnan');   % mean of force mode 3 across 5 episodes in 1 trial
        F_mode3_std=std(F_mode3,0,2,"omitnan");     % std of force mode 3 across 5 episodes in 1 trial

        F_mode4_mean = mean(F_mode4,2,'omitnan');   % mean of force mode 4 across 5 episodes in 1 trial
        F_mode4_std=std(F_mode4,0,2,"omitnan");     % std of force mode 4 across 5 episodes in 1 trial
        
% Behaviour plots of each trial averaged across 5 drift trials  

        figure(60);     % plot of total force across 5 episodes of drift
        subplot(2,5,k);
        hold on;
        x = 1:length(F_total_drift_mat);
        plot(x, F_total_drift_mat, 'LineWidth', 0.5);
        plot(x, F_total_drift_mean, 'r-', 'LineWidth', 3);
        plot(x, (F_total_drift_mean-F_total_drift_std), 'k:', 'LineWidth', 2);
        plot(x, (F_total_drift_mean+F_total_drift_std), 'k:', 'LineWidth', 2);
        xpatch = [x, fliplr(x)]; % x values are the same for both lines
        ypatch = [F_total_drift_mean-F_total_drift_std ; flip(F_total_drift_mean+F_total_drift_std)]; % y values are different for the two lines
        % Create the patch
        
        xlabel('Time');
        ylabel('% MVC');
        if contains(match_F_sensor{k}.name,"Lt","IgnoreCase",true)
            patch(xpatch, ypatch, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            title(sprintf('Total force (Left)', F_total_drift_mean, F_total_drift_std));
        else
            patch(xpatch, ypatch, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            title(sprintf('Total force (Right)', F_total_drift_mean, F_total_drift_std));  
        end

% plots for individual finger modes
        figure;             
        subplot (2,2,1);    % plot for finger mode 1
        hold on;
        x = 1:length(F_mode1);
        plot(x, F_mode1, 'LineWidth', 0.5);
        plot(x, F_mode1_mean, 'r-', 'LineWidth', 3);
        plot(x, (F_mode1_mean-F_mode1_std), 'k:', 'LineWidth', 2);
        plot(x, (F_mode1_mean+F_mode1_std), 'k:', 'LineWidth', 2);
        xpatch = [x, fliplr(x)]; % x values are the same for both lines
        ypatch = [F_mode1_mean-F_mode1_std ; flip(F_mode1_mean+F_mode1_std)]; % y values are different for the two lines
        % Create the patch
        patch(xpatch, ypatch, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        xlabel('Time');
        ylabel('%MVC');
        title(sprintf('Force mode1', F_mode1_mean, F_mode1_std));
        legend('Data', 'Mean', 'Standard Deviation');
        subplot (2,2,2);    % plot for finger mode 2
        hold on;
        x = 1:length(F_mode2);
        plot(x, F_mode2, 'LineWidth', 0.5);
        plot(x, F_mode2_mean, 'r-', 'LineWidth', 3);
        plot(x, (F_mode2_mean-F_mode2_std), 'k:', 'LineWidth', 2);
        plot(x, (F_mode2_mean+F_mode2_std), 'k:', 'LineWidth', 2);
        xpatch = [x, fliplr(x)]; % x values are the same for both lines
        ypatch = [F_mode2_mean-F_mode2_std ; flip(F_mode2_mean+F_mode2_std)]; % y values are different for the two lines
        % Create the patch
        patch(xpatch, ypatch, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        xlabel('Time');
        ylabel('%MVC');
        title(sprintf('Force mode2', F_mode2_mean, F_mode2_std));
        legend('Data', 'Mean', 'Standard Deviation');
        subplot (2,2,3);    % plot for finger mode 3
        hold on;
        x = 1:length(F_mode3);
        plot(x, F_mode3, 'LineWidth', 0.5);
        plot(x, F_mode3_mean, 'r-', 'LineWidth', 3);
        plot(x, (F_mode3_mean-F_mode3_std), 'k:', 'LineWidth', 2);
        plot(x, (F_mode3_mean+F_mode3_std), 'k:', 'LineWidth', 2);
        xpatch = [x, fliplr(x)]; % x values are the same for both lines
        ypatch = [F_mode3_mean-F_mode3_std ; flip(F_mode3_mean+F_mode3_std)]; % y values are different for the two lines
        % Create the patch
        patch(xpatch, ypatch, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        xlabel('Time');
        ylabel('%MVC');
        title(sprintf('Force mode3', F_mode3_mean, F_mode3_std));
        legend('Data', 'Mean', 'Standard Deviation');
        subplot (2,2,4);    % plot for finger mode 4
        hold on;
        x = 1:length(F_mode4);
        plot(x, F_mode4, 'LineWidth', 0.5);
        plot(x, F_mode4_mean, 'r-', 'LineWidth', 3);
        plot(x, (F_mode4_mean-F_mode4_std), 'k:', 'LineWidth', 2);
        plot(x, (F_mode4_mean+F_mode4_std), 'k:', 'LineWidth', 2);
        xpatch = [x, fliplr(x)]; % x values are the same for both lines
        ypatch = [F_mode4_mean-F_mode4_std ; flip(F_mode4_mean+F_mode4_std)]; % y values are different for the two lines
        % Create the patch
        patch(xpatch, ypatch, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        xlabel('Time');
        ylabel('%MVC');
        title(sprintf('Force mode4', F_mode4_mean, F_mode4_std));
        legend('Data', 'Mean', 'Standard Deviation');

%mode 1 & 2 plots for sensor 1
        mode_sensor1_drift_mat = cell2mat(mode_sensor1_drift(k,:));
        mode1_sensor1_drift_mat = mode_sensor1_drift_mat(:, mod(1:size(mode_sensor1_drift_mat,2),2) == 1);
        mode2_sensor1_drift_mat = mode_sensor1_drift_mat(:, mod(1:size(mode_sensor1_drift_mat,2),2) == 0);
        mode1_sensor1_drift_mean = mean(mode1_sensor1_drift_mat,2); 
        mode2_sensor1_drift_mean = mean(mode2_sensor1_drift_mat,2);
        mode1_sensor1_drift_std=std(mode1_sensor1_drift_mat,0,2);
        mode2_sensor1_drift_std=std(mode2_sensor1_drift_mat,0,2);

        figure;         % plot for mode 1 of sensor 1
        subplot(2,1,1);
        hold on;
        x = 1:length(mode1_sensor1_drift_mat);
        plot(x, mode1_sensor1_drift_mat, 'LineWidth', 0.5);
        plot(x, mode1_sensor1_drift_mean, 'r-', 'LineWidth', 3);
        plot(x, (mode1_sensor1_drift_mean-mode1_sensor1_drift_std), 'k:', 'LineWidth', 2);
        plot(x, (mode1_sensor1_drift_mean+mode1_sensor1_drift_std), 'k:', 'LineWidth', 2);
        xpatch = [x, fliplr(x)]; % x values are the same for both lines
        ypatch = [mode1_sensor1_drift_mean-mode1_sensor1_drift_std ; flip(mode1_sensor1_drift_mean+mode1_sensor1_drift_std)]; % y values are different for the two lines
        % Create the patch
        patch(xpatch, ypatch, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        xlabel('Time');
        ylabel('Mode-mag');
        title(sprintf('mode1sensor1', mode1_sensor1_drift_mean, mode1_sensor1_drift_std));
        legend('Data', 'Mean', 'Standard Deviation');

        subplot(2,1,2); % plot for mode 2 of sensor 1
        hold on;
        x = 1:length(mode2_sensor1_drift_mat);
        plot(x, mode2_sensor1_drift_mat,'LineWidth', 0.5);
        plot(x, mode2_sensor1_drift_mean, 'r-', 'LineWidth', 3);
        plot(x, (mode2_sensor1_drift_mean-mode2_sensor1_drift_std), 'k:', 'LineWidth', 2);
        plot(x, (mode2_sensor1_drift_mean+mode2_sensor1_drift_std), 'k:', 'LineWidth', 2);
        xpatch = [x, fliplr(x)]; % x values are the same for both lines
        ypatch = [mode2_sensor1_drift_mean-mode2_sensor1_drift_std ; flip(mode2_sensor1_drift_mean+mode2_sensor1_drift_std)]; % y values are different for the two lines
        % Create the patch
        patch(xpatch, ypatch, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        xlabel('Time');
        ylabel('Mode-mag');
        title(sprintf('mode2sensor1', mode2_sensor1_drift_mean, mode2_sensor1_drift_std));
        legend('Data', 'Mean', 'Standard Deviation');

%mode 1 & 2 plots for sensor 2
        mode_sensor2_drift_mat = cell2mat(mode_sensor2_drift(k,:));
        mode1_sensor2_drift_mat = mode_sensor2_drift_mat(:, mod(1:size(mode_sensor2_drift_mat,2),2) == 1);
        mode2_sensor2_drift_mat = mode_sensor2_drift_mat(:, mod(1:size(mode_sensor2_drift_mat,2),2) == 0);
        mode1_sensor2_drift_mean = mean(mode1_sensor2_drift_mat,2); 
        mode2_sensor2_drift_mean = mean(mode2_sensor2_drift_mat,2);
        mode1_sensor2_drift_std=std(mode1_sensor2_drift_mat,0,2);
        mode2_sensor2_drift_std=std(mode2_sensor2_drift_mat,0,2);


        figure;         % plot for mode 1 of sensor 2
        subplot(2,1,1);
        hold on;
        x = 1:length(mode1_sensor2_drift_mat);
        plot(x, mode1_sensor2_drift_mat, 'LineWidth', 0.5);
        plot(x, mode1_sensor2_drift_mean, 'r-', 'LineWidth', 3);
        plot(x, (mode1_sensor2_drift_mean-mode1_sensor2_drift_std), 'k:', 'LineWidth', 2);
        plot(x, (mode1_sensor2_drift_mean+mode1_sensor2_drift_std), 'k:', 'LineWidth', 2);
        xpatch = [x, fliplr(x)]; % x values are the same for both lines
        ypatch = [mode1_sensor2_drift_mean-mode1_sensor2_drift_std ; flip(mode1_sensor2_drift_mean+mode1_sensor2_drift_std)]; % y values are different for the two lines
        % Create the patch
        patch(xpatch, ypatch, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        xlabel('Time');
        ylabel('Mode-mag');
        title(sprintf('mode1sensor2', mode1_sensor2_drift_mean, mode1_sensor2_drift_std));
        legend('Data', 'Mean', 'Standard Deviation');

        subplot(2,1,2); % plot for mode 2 of sensor 2
        hold on;
        x = 1:length(mode2_sensor2_drift_mat);
        plot(x, mode2_sensor2_drift_mat,'LineWidth', 0.5);
        plot(x, mode2_sensor2_drift_mean, 'r-', 'LineWidth', 3);
        plot(x, (mode2_sensor2_drift_mean-mode2_sensor2_drift_std), 'k:', 'LineWidth', 2);
        plot(x, (mode2_sensor2_drift_mean+mode2_sensor2_drift_std), 'k:', 'LineWidth', 2);
        xpatch = [x, fliplr(x)]; % x values are the same for both lines
        ypatch = [mode2_sensor2_drift_mean-mode2_sensor2_drift_std ; flip(mode2_sensor2_drift_mean+mode2_sensor2_drift_std)]; % y values are different for the two lines
        % Create the patch
        patch(xpatch, ypatch, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        xlabel('Time');
        ylabel('Mode-mag');
        title(sprintf('mode2sensor2', mode2_sensor2_drift_mean, mode2_sensor2_drift_std));
        legend('Data', 'Mean', 'Standard Deviation');

%mode 1 & 2 plots for sensor 12
        mode_sensor12_drift_mat = cell2mat(mode_sensor12_drift(k,:));
        mode1_sensor12_drift_mat = mode_sensor12_drift_mat(:, mod(1:size(mode_sensor12_drift_mat,2),2) == 1);
        mode2_sensor12_drift_mat = mode_sensor12_drift_mat(:, mod(1:size(mode_sensor12_drift_mat,2),2) == 0);
        mode1_sensor12_drift_mean = mean(mode1_sensor12_drift_mat,2); 
        mode2_sensor12_drift_mean = mean(mode2_sensor12_drift_mat,2);
        mode1_sensor12_drift_std=std(mode1_sensor12_drift_mat,0,2);
        mode2_sensor12_drift_std=std(mode2_sensor12_drift_mat,0,2);

        figure;         % plot for mode 1 of sensor 12
        subplot(2,1,1);
        hold on;
        x = 1:length(mode1_sensor12_drift_mat);
        plot(x, mode1_sensor12_drift_mat, 'LineWidth', 0.5);
        plot(x, mode1_sensor12_drift_mean, 'r-', 'LineWidth', 3);
        plot(x, (mode1_sensor12_drift_mean-mode1_sensor12_drift_std), 'k:', 'LineWidth', 2);
        plot(x, (mode1_sensor12_drift_mean+mode1_sensor12_drift_std), 'k:', 'LineWidth', 2);
        xpatch = [x, fliplr(x)]; % x values are the same for both lines
        ypatch = [mode1_sensor12_drift_mean-mode1_sensor12_drift_std ; flip(mode1_sensor12_drift_mean+mode1_sensor12_drift_std)]; % y values are different for the two lines
        % Create the patch
        patch(xpatch, ypatch, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        xlabel('Time');
        ylabel('Mode-mag');
        title(sprintf('mode1sensor12', mode1_sensor12_drift_mean, mode1_sensor12_drift_std));
        legend('Data', 'Mean', 'Standard Deviation');

        subplot(2,1,2); % plot for mode 2 of sensor 12
        hold on;
        x = 1:length(mode2_sensor12_drift_mat);
        plot(x, mode2_sensor12_drift_mat,'LineWidth', 0.5);
        plot(x, mode2_sensor12_drift_mean, 'r-', 'LineWidth', 3);
        plot(x, (mode2_sensor12_drift_mean-mode2_sensor12_drift_std), 'k:', 'LineWidth', 2);
        plot(x, (mode2_sensor12_drift_mean+mode2_sensor12_drift_std), 'k:', 'LineWidth', 2);
        xpatch = [x, fliplr(x)]; % x values are the same for both lines
        ypatch = [mode2_sensor12_drift_mean-mode2_sensor12_drift_std ; flip(mode2_sensor12_drift_mean+mode2_sensor12_drift_std)]; % y values are different for the two lines
        % Create the patch
        patch(xpatch, ypatch, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        xlabel('Time');
        ylabel('Mode-mag');
        title(sprintf('mode2sensor12', mode2_sensor12_drift_mean, mode2_sensor12_drift_std));
        legend('Data', 'Mean', 'Standard Deviation');
        
end

%% UCM analysis
 for t=1:size(F_drift,1) % analysing each trial

     % UCM analysis of difference between post and pre for sensor1
     d1_s1=d_M1_sensor1(:,:,t);
     d2_s1=d_M2_sensor1(:,:,t);
     d_s1=[d1_s1(:),d2_s1(:)];
     ucm_d_s1(t,:)=UCM(d_s1,jacob_sensor1(t,:));
     display(ucm_d_s1(t,:));

     % UCM analysis of difference between post and pre for sensor2
     d1_s2=d_M1_sensor2(:,:,t);
     d2_s2=d_M2_sensor2(:,:,t);
     d_s2=[d1_s2(:),d2_s2(:)];
     ucm_d_s2(t,:)=UCM(d_s2,jacob_sensor2(t,:));
     display(ucm_d_s2(t,:));

     % UCM analysis of difference between post and pre for sensor1+2
     d1_s12=d_M1_sensor12(:,:,t);
     d2_s12=d_M2_sensor12(:,:,t);
     d_s12=[d1_s12(:),d2_s12(:)];
     ucm_d_s12(t,:)=UCM(d_s12,jacob_sensor12(t,:));
     display(ucm_d_s12(t,:));

     % UCM analysis of pre for sensor1
     pre1_s1=pre_M1_sensor1(:,:,t);
     pre2_s1=pre_M2_sensor1(:,:,t);
     pre_s1=[pre1_s1(:),pre2_s1(:)];
     ucm_pre_s1(t,:)=UCM(pre_s1,jacob_sensor1(t,:));
     display(ucm_pre_s1(t,:));

     % UCM analysis of pre for sensor2
     pre1_s2=pre_M1_sensor2(:,:,t);
     pre2_s2=pre_M2_sensor2(:,:,t);
     pre_s2=[pre1_s2(:),pre2_s2(:)];
     ucm_pre_s2(t,:)=UCM(pre_s2,jacob_sensor2(t,:));
     display(ucm_pre_s2(t,:));

     % UCM analysis of pre for sensor1+2
     pre1_s12=pre_M1_sensor12(:,:,t);
     pre2_s12=pre_M2_sensor12(:,:,t);
     pre_s12=[pre1_s12(:),pre2_s12(:)];
     ucm_pre_s12(t,:)=UCM(pre_s12,jacob_sensor12(t,:));
     display(ucm_pre_s12(t,:));

     % UCM analysis of post for sensor1
     post1_s1=post_M1_sensor1(:,:,t);
     post2_s1=post_M2_sensor1(:,:,t);
     post_s1=[post1_s1(:),post2_s1(:)];
     ucm_post_s1(t,:)=UCM(post_s1,jacob_sensor1(t,:));
     display(ucm_post_s1(t,:));

     % UCM analysis of post for sensor2
     post1_s2=post_M1_sensor2(:,:,t);
     post2_s2=post_M2_sensor2(:,:,t);
     post_s2=[post1_s2(:),post2_s2(:)];
     ucm_post_s2(t,:)=UCM(post_s2,jacob_sensor2(t,:));
     display(ucm_post_s2(t,:));

     % UCM analysis of post for sensor1+2
     post1_s12=post_M1_sensor12(:,:,t);
     post2_s12=post_M2_sensor12(:,:,t);
     post_s12=[post1_s12(:),post2_s12(:)];
     ucm_post_s12(t,:)=UCM(post_s12,jacob_sensor12(t,:));
     display(ucm_post_s12(t,:));

 end

 %% Synergy analysis for Force mode
 pre_F_mode_lt=nan((L_match*n_drift),4);
 post_F_mode_lt=nan((L_match*n_drift),4);
 pre_F_mode_rt=nan((L_match*n_drift),4);
 post_F_mode_rt=nan((L_match*n_drift),4);
z=1;
  for i=1:size(F_drift,1)
     
     if contains(match_F_sensor{i}.name,"Lt","IgnoreCase",true)
         pre_F_mode_lt(((i-1)*5+1):(i*5),:) = pre_F_mode(:,:,i);
         post_F_mode_lt(((i-1)*5+1):(i*5),:) = post_F_mode(:,:,i);

     else
         
         pre_F_mode_rt(((z-1)*5+1):(z*5),:) = pre_F_mode(:,:,i);
         post_F_mode_rt(((z-1)*5+1):(z*5),:) = post_F_mode(:,:,i);

         z=z+1;
     end
     
 end

 ucm_pre_F_mode_lt=UCM(pre_F_mode_lt,ones(1,4));
 ucm_pre_F_mode_rt=UCM(pre_F_mode_rt,ones(1,4));

 ucm_post_F_mode_lt=UCM(post_F_mode_lt,ones(1,4));
 ucm_post_F_mode_rt=UCM(post_F_mode_rt,ones(1,4));