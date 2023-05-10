%% Read all the force data files
clc
clear
close all
subjid='subj08';

q=cd("Force_drift");
p=cd(subjid);
FileList=dir('*MVC.txt');
L=length(FileList);
data=cell(L,1);
for i=1:L
data{i}.name=FileList(i).name;

finger_forces=load(FileList(i).name);

data{i}.signals=max(finger_forces);
end

cd(fullfile(q));
mkdir(subjid);
cd(subjid);
save('MVC_data.mat','data');
cd(fullfile(p,subjid));



%%

FileList=dir('*enslaving*.txt');
L=length(FileList);
data=cell(L,1);
for i=1:L
data{i}.name=FileList(i).name;

finger_forces=load(FileList(i).name);

data{i}.signals=finger_forces(:,:);

end

cd(fullfile(q));
mkdir(subjid);
cd(subjid);
save('enslaving_data.mat','data')
cd(fullfile(p,subjid));


%%

FileList=dir('*IMRL_raw*.txt');
L=length(FileList);
data=cell(L,1);
for i=1:L
data{i}.name=FileList(i).name;

finger_forces=load(FileList(i).name);

data{i}.signals=finger_forces(:,:);

end

cd(fullfile(q));
mkdir(subjid);
cd(subjid);
save('drift_data.mat','data')
cd(fullfile(q));

