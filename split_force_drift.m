%% Split function for force data where all the episodes are individually cut and each individual cycles are analyzed for acceptance in the cyclical episode
function [F_drift,T_drift,FT_drift]=split_force_drift(F_data,filename,ENSL_L,ENSL_R,MVC)

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



ind_cycle_start=find(F_data(:,1)<5,1,"last");    % start index of 1st period (cyclical)
ind_cycle_end=find(F_data(:,1)<20000,1,"last");  % end index of 1st period (cyclical)



figure

subplot(3,1,1)
plot(F_data(:,1),sum(F_data(:,2:end),2),'DisplayName','data(:,2:end)'); % plotting the force behaviour for the whole trial
xline([28000 51000 74000 97000 120000],"Color",'g');
xline([26750 49750 72750 95750 118750],"Color",'r');

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




%% chopping each period of drift in 1 trial
k=26750
l=43000
    for j=1:5
        F_drift{j}=F_data(find(F_data(:,1)>k,1,"first"):find(F_data(:,1)<l,1,"last"),2:end); % store force data from variable "data" for all 4 force but from the peak location identefied to the next peak location
        T_drift{j}=F_data(find(F_data(:,1)>k,1,"first"):find(F_data(:,1)<l,1,"last"),1);
        FT_drift{j}=F_data(find(F_data(:,1)>k,1,"first"):find(F_data(:,1)<l,1,"last"),:);

        k=k+23000;
        l=l+23000;
    end
