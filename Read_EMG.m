%% EMG data
clc
clear
close all
subjid='subj08';


q=cd("EMG_drift");
p=cd(subjid);
r=cd('txt');
D1 = dir("*Sensor 1*") ;
D2 = dir("*Sensor 2*") ;
L1=length(D1); 
L2=length(D2); 
sensor1_data=cell(L1,1);
sensor2_data=cell(L2,1);

%% sensor1 = Flexor
for k=1:length(D1) 
currentD=D1(k).name;
cd(currentD)
FileList=dir('*MFR.txt');
sensor1_data{k}.name=FileList.name;
sensor1_data{k}.signals=importdata(FileList.name). data;
cd(fullfile(r,'txt'));
end


%% sensor2 = Extensor

for k=1:length(D2)
currentD=D2(k).name;
cd(currentD)
FileList=dir('*MFR.txt');
sensor2_data{k}.name=FileList.name;
sensor2_data{k}.signals=importdata(FileList.name).data;
cd(fullfile(r,'txt'));
end

%% 
cd(fullfile(q));
mkdir(subjid);
cd(subjid);
save('sensor1.mat','sensor1_data');
save('sensor2.mat','sensor2_data');
cd(fullfile(q));