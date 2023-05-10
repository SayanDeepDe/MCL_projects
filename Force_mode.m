clc
clear
close all
subjid='subj10';
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



%% enslaving coefficient
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

%% split = splitting the episodes
t_drift=cell(length(drift_f),1);
f_drift=cell(length(drift_f),1);
for i=1:length(t_drift)
    t_drift{i}.signals=drift_f{i}.signals(:,1);
    t_drift{i}.name=drift_f{i}.name;
    f_drift{i}.signals=drift_f{i}.signals(:,2:5);
    f_drift{i}.name=drift_f{i}.name;
end

%% Entering each trial and processing them individually

T_drift=cell(10,5);
F_drift=cell(10,5);
FT_drift=cell(10,5);



for j=1:length(t_drift) % each trial
   [F_d,T_d,FT_d]=split_force_drift(drift_f{j}.signals,drift_f{j}.name,ENSL_L,ENSL_R,MVC) ;


%   storing chopped drift trials. Each row=trial, column=episode
    [F_drift{j,:}]=F_d{:};                  %storing Force episodes from each trial
    [T_drift{j,:}]=T_d{:};                  %storing Time(Force) episodes from each trial

end


%% converting Force data in drift episodes to Force-mode


 for i=1:length(F_drift)
     
     if contains(drift_f{i}.name,"Lt","IgnoreCase",true)
         for j=1:size(F_drift,2)
                 [F_datamode{i,j}]=deal((inv(ENSL_L)* (F_drift{i,j})')');
         end
     else
         for j=1:size(F_drift,2)
                [F_datamode{i,j}]=deal((inv(ENSL_R)* (F_drift{i,j})')'); 
         end
     end
     
 end


%% Behaviour plot

pre_F_mode=nan(5,4,10);
post_F_mode=nan(5,4,10);

for k=1:length(F_drift) % moving through each trial.
    pre_F=26750;    % start of drift episode 1
    post_F=43000;   % end of drift episode 1

    for l=1:size(F_drift,2) % moving through each episode of drift of a single trial.
        F_total_drift{k,l}=sum(F_drift{k,l},2);
                % start indexes of 3 time windows pre drift
        pre_F_1=pre_F;
        post_F_1=post_F-1250;

       for w=1:3 % moving through each time window
           [~,pre_T_F(1,w)]=min(abs(T_drift{k,l}-pre_F_1)); 
           [~,pre_T_F(2,w)]=min(abs(T_drift{k,l}-(pre_F_1+250)));

           [~,post_T_F(1,w)]=min(abs(T_drift{k,l}-(post_F_1)));
           [~,post_T_F(2,w)]=min(abs(T_drift{k,l}-(post_F_1+250)));

           pre_F_1=pre_F_1+500;
           post_F_1=post_F_1+500;
       end

               %Forcemode
           pre_F_mode(l,:,k)=mean((F_datamode{k,l}(pre_T_F(1,3):pre_T_F(2,3),:)));
           post_F_mode(l,:,k)=mean((F_datamode{k,l}(post_T_F(1,3):post_T_F(2,3),:)));
    end

     % Horizontally concatenate, pad with NaNs
        maxNumRows = max(cellfun(@(c) length(c), F_total_drift(k,:)));  % max number of columns
        F_total_drift_mat = cell2mat(cellfun(@(c){padarray(c,[maxNumRows-length(c),0],NaN,'Post')}, F_total_drift(k,:)));
        T_drift_mat = cell2mat(cellfun(@(c){padarray(c,[maxNumRows-length(c),0],NaN,'Post')}, T_drift(k,:)));
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
        if (k<6)
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


end

%% Synergy analysis for Force mode
z=1;

  for i=1:length(F_drift)
     
     if contains(drift_f{i}.name,"Lt","IgnoreCase",true)
         pre_F_mode_lt(((i-1)*5+1):(i*5),:) = pre_F_mode(:,:,i);
         post_F_mode_lt(((i-1)*5+1):(i*5),:) = post_F_mode(:,:,i);

     elseif contains(drift_f{i}.name,"Rt","IgnoreCase",true)
         
         pre_F_mode_rt(((z-1)*5+1):(z*5),:) = pre_F_mode(:,:,i);
         post_F_mode_rt(((z-1)*5+1):(z*5),:) = post_F_mode(:,:,i);

         z=z+1;
     end
     
 end

 ucm_pre_F_mode_lt=UCM(pre_F_mode_lt,ones(1,4));
 ucm_pre_F_mode_rt=UCM(pre_F_mode_rt,ones(1,4));

 ucm_post_F_mode_lt=UCM(post_F_mode_lt,ones(1,4));
 ucm_post_F_mode_rt=UCM(post_F_mode_rt,ones(1,4));