%% This functio is used to split the force, EMG files into episodes of cycles and drifts.
% It cuts individual cycles, applies the criterion to accept the cyles and
% make 20 phases for each cycle.
function [sensor1_cycle,sensor2_cycle,T_sensor_cycle,sensor1_cycle_ave,sensor1_cycle_ave_new,sensor2_cycle_ave,sensor2_cycle_ave_new,T_cycle_ave,T_sensor_cycle_ave,F_cycle_ave,F_cycle_tot_ave,F_cycle_tot_ave_new,T_sensor_drift,sensor1_drift,sensor2_drift,F_drift,T_drift,FT_drift,sensor1_T_drift,sensor2_T_drift]=split_drift(sensor1,sensor2,F_data,filename,ENSL_L,ENSL_R,MVC)

Lt_IMRL=contains(filename,"Lt_IMRL","IgnoreCase",true);
Rt_IMRL=contains(filename,"Rt_IMRL","IgnoreCase",true);

[b,a] = butter(2,10/(200/2));
F_data(:,2:5)=filter(b,a,F_data(:,2:5)); % filtering the data = Butterworth filter
%% Normalizing the data wrt MVC
for i=1:length(MVC)
    if Lt_IMRL==1 && contains(MVC{i}.name,"Lt_IMRL","IgnoreCase",true)==1
    F_data(:,2:5)=F_data(:,2:5)/MVC{i}.signals;
    elseif Rt_IMRL==1 && contains(MVC{i}.name,"Rt_IMRL","IgnoreCase",true)==1
    F_data(:,2:5)=F_data(:,2:5)/MVC{i}.signals;
    end
end


%% 
ind_cycle_start=find(F_data(:,1)<5,1,"last");    % start index of 1st period (cyclical)
ind_cycle_end=find(F_data(:,1)<20000,1,"last");  % end index of 1st period (cyclical)

%% Behavorial plot
figure  % behavorial plot for entire trial with the time phases for pre and post drift analysis

subplot(3,1,1)
plot(F_data(:,1),sum(F_data(:,2:end),2),'DisplayName','data(:,2:end)'); % plotting the force behaviour for the whole trial
xline([28000 51000 74000 97000 120000],"Color",'g');
xline([26750 49750 72750 95750 118750],"Color",'r');

subplot(3,1,2)
plot(sensor1(:,1),sensor1(:,2:end),'DisplayName','data(:,2:end)'); % plotting the EMG-flexor behaviour for the whole trial

subplot(3,1,3)
plot(sensor2(:,1),sensor2(:,2:end),'DisplayName','data(:,2:end)'); % plotting the EMG-extensor behaviour for the whole trial
%% Finding peak in half cycles without up and down
t_data=F_data(ind_cycle_start:ind_cycle_end,1); % time of the 1st period (cyclical)
data_cycle=sum(F_data(ind_cycle_start:ind_cycle_end,2:5),2);
data_mean=data_cycle-mean(data_cycle);
[~,ind_peak,width_mean,height_mean] = findpeaks(data_mean); % find all peak and the amplitude of the peaks from the data
f_peak=data_cycle(ind_peak);% finding the force values at the peak locations
ind_true_peak=ind_peak(width_mean>100);
f_true_peak=f_peak(width_mean>100);
% ind_true_peak=ind_peak(height_mean>15); % only considering the peak values that has an amplitude of greater than 1
% f_true_peak=f_peak(height_mean>15); % also considering the force values at only amplitude > 1
% ind_true_peak=ind_peak(height_mean>15);
% f_true_peak=f_peak(height_mean>15);

figure

subplot(5,1,1)
plot(t_data,data_cycle)
hold on
plot(t_data,F_data(ind_cycle_start:ind_cycle_end,2:5))
hold on
plot(t_data(ind_true_peak),f_true_peak,'*')
title('Total Force & individual finger force')
subtitle('Peaks on total force')

%% chopping each cycle force and respective EMG data

%  make 2 matrix containing the Force data and the EMG data for the cycle 
cyc=1; % initialising the index for 1st cycle
for i=1:length(ind_true_peak)-1

    F_cycle{cyc}=F_data(ind_true_peak(i):ind_true_peak(i+1),2:5); % store force data from variable "data" for all 4 force but from the peak location identefied to the next peak location
    T_cycle{cyc}=F_data(ind_true_peak(i):ind_true_peak(i+1),1);
    FT_cycle{cyc}=F_data(ind_true_peak(i):ind_true_peak(i+1),:);

    %translating index from time of force matrix to time of EMG matrix
    s=find(sensor1(:,1)<round(T_cycle{cyc}(1,1)/1000,2),1,'last');                  % start index of EMG data
    e=find(sensor1(:,1)<round(T_cycle{cyc}(size(T_cycle{cyc},1))/1000,2),1,'last'); % end index of EMG data

    sensor1_cycle{cyc}=sensor1(s:e,2:size(sensor1,2));
    sensor1_T_cycle{cyc}=sensor1(s:e,:);

    sensor2_cycle{cyc}=sensor2(s:e,2:size(sensor2,2));
    sensor2_T_cycle{cyc}=sensor2(s:e,:);

    T_sensor_cycle{cyc}=sensor1(s:e,1); %the x-axis from EMG data

    cyc=cyc+1; % to go to the next location of the matrix and store data for next cycle
end

%% omitting the unaccepted cycles

for i=1:size(F_cycle,2)
    l(i)=length(F_cycle{i}); % taking the length/width of each cycle
end

d_l=find(l<(mean(l)-0.3*mean(l))|l>(mean(l)+0.3*mean(l))); % ommiting cycles whose length is > or < than 30% of mean length of each cycle
F_cycle(d_l)=[];
T_cycle(d_l)=[];
FT_cycle(d_l)=[];
sensor1_cycle(d_l)=[];
sensor1_T_cycle(d_l)=[];
sensor2_cycle(d_l)=[];
sensor2_T_cycle(d_l)=[];
T_sensor_cycle(d_l)=[];

%% slicing the force and EMG cycles into 50 equal parts
n_bins=20; % initialising the number of bins to be created

for j=1:size(F_cycle,2)
    F_bins=array_split(F_cycle{j},n_bins);
    T_bins=array_split(T_cycle{j},n_bins);
    sensor1_bins=array_split(sensor1_cycle{j},n_bins);
    sensor2_bins=array_split(sensor2_cycle{j},n_bins);
    T_sensor_bins=array_split(T_sensor_cycle{j},n_bins);

    F_cycle_bins=nan(n_bins,(size(F_data,2)-1));
    T_cycle_bins=nan(n_bins,1);
    sensor1_cycle_bins=nan(n_bins,(size(sensor1,2)-1));
    sensor2_cycle_bins=nan(n_bins,(size(sensor2,2)-1));
    T_sensor_cycle_bins=nan(n_bins,1);

    for k=1:size(F_bins) % average across time windows
        F_cycle_bins(k,:)=mean((F_bins{k}),1);
        T_cycle_bins(k,:)=mean((T_bins{k}),1);
        sensor1_cycle_bins(k,:)=mean((sensor1_bins{k}),1);
        sensor2_cycle_bins(k,:)=mean((sensor2_bins{k}),1);
        T_sensor_cycle_bins(k,:)=mean((T_sensor_bins{k}),1);
    end

    F_cycle_ave{j}=F_cycle_bins;
    T_cycle_ave{j}=T_cycle_bins;
    sensor1_cycle_ave{j}=sensor1_cycle_bins;
    sensor2_cycle_ave{j}=sensor2_cycle_bins;
    T_sensor_cycle_ave{j}=T_sensor_cycle_bins;

end

%% Average across 20-phased cycle for NEW jacobian analysis
sensor1_cycle_ave_mat=cat(3,sensor1_cycle_ave{:});
sensor1_cycle_ave_new=mean(sensor1_cycle_ave_mat,3);
sensor2_cycle_ave_mat=cat(3,sensor2_cycle_ave{:});
sensor2_cycle_ave_new=mean(sensor2_cycle_ave_mat,3);
F_cycle_tot_ave_across=cellfun(@(x) mean(x, 2), F_cycle_ave, 'UniformOutput', false);
F_cycle_tot_ave_mat=cat(3,F_cycle_tot_ave_across{:});
F_cycle_tot_ave_new=mean(F_cycle_tot_ave_mat,3);
%% concateting the accepted cycles together
F_cycle=cat(1,F_cycle{:}); 
F_cycle_tot=sum(F_cycle,2);
T_cycle=cat(1,T_cycle{:});
sensor1_cycle=cat(1,sensor1_cycle{:}); 
sensor2_cycle=cat(1,sensor2_cycle{:});
T_sensor_cycle=cat(1,T_sensor_cycle{:});

F_cycle_ave=cat(1,F_cycle_ave{:}); 
F_cycle_tot_ave=sum(F_cycle_ave,2);
T_cycle_ave=cat(1,T_cycle_ave{:});
sensor1_cycle_ave=cat(1,sensor1_cycle_ave{:}); 
sensor2_cycle_ave=cat(1,sensor2_cycle_ave{:});
T_sensor_cycle_ave=cat(1,T_sensor_cycle_ave{:});

%% Finger mode computation for cyclical episode
if contains(filename,"Lt","IgnoreCase",true)
        M_data_cycle=(inv(ENSL_L) * (F_cycle)')';
   
      
     else
        M_data_cycle=(inv(ENSL_R) * (F_cycle)')'; 
       
end

%% plotting various data into 1 figure
subplot(5,1,2)
plot(T_cycle,F_cycle_tot)
hold on
plot(t_data(ind_true_peak),f_true_peak,'*')
hold on 
plot(T_cycle,F_cycle)
title('Total Force & individual finger force')
subtitle('After excluding the unaccepted cycles')


subplot(5,1,3)
plot(T_cycle,sum(M_data_cycle,2))
hold on
plot(T_cycle,M_data_cycle)
hold on
plot(t_data(ind_true_peak),f_true_peak,'*')
title('Total Mode Force & individual finger-mode force')
subtitle('After excluding the unaccepted cycles')

subplot(5,1,4)
plot(T_sensor_cycle,sensor1_cycle)
title('Sensor1 MU = Flexor')
subtitle('Aligned with the accepted force cycles')

subplot(5,1,5)
plot(T_sensor_cycle,sensor2_cycle)
title('Sensor2 MU = Extensor')
subtitle('Aligned with the accepted force cycles')

%% chopping each period of drift in 1 trial
k=26750;
l=43000;
    for j=1:5
        F_drift{j}=F_data(find(F_data(:,1)>k,1,"first"):find(F_data(:,1)<l,1,"last"),2:end); % store force data from variable "data" for all 4 force but from the peak location identefied to the next peak location
        T_drift{j}=F_data(find(F_data(:,1)>k,1,"first"):find(F_data(:,1)<l,1,"last"),1);
        FT_drift{j}=F_data(find(F_data(:,1)>k,1,"first"):find(F_data(:,1)<l,1,"last"),:);
        sensor1_drift{j}=sensor1(find(sensor1(:,1)>(k/1000),1,"first"):find(sensor1(:,1)<(l/1000),1,"last"),2:end);
        sensor1_T_drift{j}=sensor1(find(sensor1(:,1)>(k/1000),1,"first"):find(sensor1(:,1)<(l/1000),1,"last"),:);
        sensor2_drift{j}=sensor2(find(sensor2(:,1)>(k/1000),1,"first"):find(sensor2(:,1)<(l/1000),1,"last"),2:end);
        sensor2_T_drift{j}=sensor2(find(sensor2(:,1)>(k/1000),1,"first"):find(sensor2(:,1)<(l/1000),1,"last"),:);
        T_sensor_drift{j}=sensor1(find(sensor1(:,1)>(k/1000),1,"first"):find(sensor1(:,1)<(l/1000),1,"last"),1);
        k=k+23000;
        l=l+23000;
    end

% T_sensor_drift=cat(1,T_sensor_drift{:});
% sensor1_drift=cat(1,sensor1_drift{:});
% sensor2_drift=cat(1,sensor2_drift{:});
% T_drift=cat(1,T_drift{:});
% F_drift=cat(1,F_drift{:});

% figure(3)
% subplot(3,1,1)
% plot(T_drift,F_drift)
% xline([28000 51000 74000 97000 120000],"Color",'g');
% xline([27000 50000 73000 96000 119000],"Color",'r');
% 
% subplot(3,1,2)
% plot(T_sensor_drift,sensor1_drift)
% xline([28 51 74 97 120],"Color",'g');
% xline([27 50 73 96 119],"Color",'r');
% 
% subplot(3,1,3)
% plot(T_sensor_drift,sensor2_drift)
% xline([28 51 74 97 120],"Color",'g');
% xline([27 50 73 96 119],"Color",'r');
     
        

%     A = repmat(T_cyc,[1 length(T_sensor_cyc)]);
%     [minValue,closestIndex] = min(abs(A-(T_sensor_cyc*1000)')); %finding the closest value and the matching index of T_sensor_cyc(EMG time) in T_cyc(Force time)
%     closestValue = T_cyc(closestIndex) ;
%     F_cyc_tot_ds = F_cyc_tot(closestIndex);
% 
%     s=find(sensor1(:,1)<round(T_cycle{cyc}(1,1)/1000,2),1,'last'); %translating index from time of force matrix to time of EMG matrix
%     e=find(sensor1(:,1)<round(T_cycle{cyc}(size(T_cycle{cyc},1))/1000,2),1,'last');
%     sensor1_cycle{cyc}=sensor1(s:e,2:size(sensor1,2));

