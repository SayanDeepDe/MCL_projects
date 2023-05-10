function [k_I,k_M,k_R,k_L,lt,rt,n]=enslave_coeff(data,filename)
%%This function takes in the data and the filename and calculate the
%%enslaving coefficient for that task finger.

timestamp=data(:,1); %taking the timestamps out from the 1st column of raw data
ind_t1_start=find(timestamp<3500,1,"last"); %trial 1 starts from 3500ms
ind_t1_end=find(timestamp<5500,1,"last"); %trial 1 ends at 5500ms
ind_t2_start=find(timestamp<15500,1,"last"); %trial 2 starts from 15500ms
ind_t2_end=find(timestamp<17500,1,"last"); %trial 2 ends at 17500ms
data_t1=data(ind_t1_start:ind_t1_end,2:5); %extracting data for trial 1 from 3500-5500ms
data_t2=data(ind_t2_start:ind_t2_end,2:5); %extracting data for trial 2 from 15500-17500ms
F_total_t1=sum(data_t1,2); %total force (IMRL) for trial 1
F_total_t2=sum(data_t2,2); %total force (IMRL) for trial 2



% adding value of 1 to detect the hand and task finger
Lt_I=contains(filename,"Lt_I","IgnoreCase",true);
Lt_M=contains(filename,"Lt_M","IgnoreCase",true);
Lt_R=contains(filename,"Lt_R","IgnoreCase",true);
Lt_L=contains(filename,"Lt_L","IgnoreCase",true);
Rt_I=contains(filename,"Rt_I","IgnoreCase",true);
Rt_M=contains(filename,"Rt_M","IgnoreCase",true);
Rt_R=contains(filename,"Rt_R","IgnoreCase",true);
Rt_L=contains(filename,"Rt_L","IgnoreCase",true);

%initialising the values to 0 which changes to 1 depending on the hand 
lt=0;
rt=0;

%% enslaving coeeff by regressing the individual finger force against total force
if (Lt_I==1 || Lt_M==1 || Lt_R==1 || Lt_L==1) % for left hand task
[k_I(1,1),~,~,~,k_I(1,2:5)]=regress(data_t1(:,1),F_total_t1); %coeff and stats in the 1st row for trial 1 for I finger
[k_I(2,1),~,~,~,k_I(2,2:5)]=regress(data_t2(:,1),F_total_t2); %coeff and stats in the 1st row for trial 2 for I finger
[k_M(1,1),~,~,~,k_M(1,2:5)]=regress(data_t1(:,2),F_total_t1); %same for M finger trial 1
[k_M(2,1),~,~,~,k_M(2,2:5)]=regress(data_t2(:,2),F_total_t2); %same for M finger trial 2
[k_R(1,1),~,~,~,k_R(1,2:5)]=regress(data_t1(:,3),F_total_t1); %same for R finger trial 1
[k_R(2,1),~,~,~,k_R(2,2:5)]=regress(data_t2(:,3),F_total_t2); %same for R finger trial 2
[k_L(1,1),~,~,~,k_L(1,2:5)]=regress(data_t1(:,4),F_total_t1); %same for L finger trial 1
[k_L(2,1),~,~,~,k_L(2,2:5)]=regress(data_t2(:,4),F_total_t2); %same for L finger trial 2
lt=1; %lt=1 when the task hand is left

else %for right hand task
[k_I(1,1),~,~,~,k_I(1,2:5)]=regress(data_t1(:,4),F_total_t1);
[k_I(2,1),~,~,~,k_I(2,2:5)]=regress(data_t2(:,4),F_total_t2);
[k_M(1,1),~,~,~,k_M(1,2:5)]=regress(data_t1(:,3),F_total_t1);
[k_M(2,1),~,~,~,k_M(2,2:5)]=regress(data_t2(:,3),F_total_t2);
[k_R(1,1),~,~,~,k_R(1,2:5)]=regress(data_t1(:,2),F_total_t1);
[k_R(2,1),~,~,~,k_R(2,2:5)]=regress(data_t2(:,2),F_total_t2);
[k_L(1,1),~,~,~,k_L(1,2:5)]=regress(data_t1(:,1),F_total_t1);
[k_L(2,1),~,~,~,k_L(2,2:5)]=regress(data_t2(:,1),F_total_t2); 
rt=1; %rt=1 when the task hand is right

end

%% comparing the R-sqr between the 2 trials of the task finger
if (Lt_I==1 || Rt_I==1) % for index finger task
    n=1; % choosing the column in the enslaving matrix based on task finger
    if k_I(1,2)>k_I(2,2) % choosing the trial depending on highest R-sqr
        trial=1;
    else
        trial=2;
    end

elseif (Lt_M==1 || Rt_M==1) % for middle finger task
    n=2; % n=1(index), 2(middle), 3(ring), 4(little)
    if k_M(1,2)>k_M(2,2)
        trial=1;
    else
        trial=2;
    end
       
elseif (Lt_R==1 || Rt_R==1) % for ring finger task
    n=3;
    if k_R(1,2)>k_R(2,2)
        trial=1;
    else
        trial=2;
    end

elseif (Lt_L==1 || Rt_L==1) % for little finger task
    n=4;
    if k_L(1,2)>k_L(2,2)
        trial=1;
    else
        trial=2;
    end    
end

%% coeff based on the higher R-sqr trial
k_I=k_I(trial,1);
k_M=k_M(trial,1);
k_R=k_R(trial,1);
k_L=k_L(trial,1);

