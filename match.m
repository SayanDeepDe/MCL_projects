function [match_F_sensor,match_sensor1,match_sensor2, null_sensor1, null_sensor2, L_match]=match(drift_f,sensor1,sensor2,L_drift,L_sensor1,L_sensor2,critical_n_MU)

%% checking empty trials for sensor1
n = 5; % replace 5 with the maximum value of i in your loop
k=0;
results_sensor1 = zeros(1, L_drift); % initialize the vector of results to all zeros

for j = 1:L_sensor1
    filename = sensor1{j}.name;
    found_i = false; % initialize a flag to indicate whether we've found i in the filename
    for i = 1:n
        if contains(filename, ['_Lt_t' num2str(i)])
            results_sensor1(i) = 1; % set the ith element of the results vector to 1
            found_i = true; % set the flag to true
            k=i;
            break; % break out of the inner loop once we've found i

        elseif contains(filename, ['_Rt_t' num2str(i)])
            results_sensor1(i+5) = 1; % set the ith element of the results vector to 1
            found_i = true; % set the flag to true
            k=i;
            break; % break out of the inner loop once we've found i

        elseif i>k && contains(filename, ['_Lt_t' num2str(i)])==0
            results_sensor1(i) = 0; % set the ith element of the results vector to 1

        elseif i>k && contains(filename, ['_Rt_t' num2str(i)])==0
            results_sensor1(i+5) = 0; % set the ith element of the results vector to 1
        end
    end
    if ~found_i
        disp(['The filename "', filename, '" does not contain any of the numbers 1 to ', num2str(n)])
    end
    if all(results_sensor1) % check if all elements of results are nonzero
        disp('All desired numbers have been found!')
        break; % break out of the outer loop if we've found all the numbers
    end
end

%% checking empty trials for sensor2
n = 5; % replace 5 with the maximum value of i in your loop
k=0;
results_sensor2 = zeros(1, L_drift); % initialize the vector of results to all zeros

for j = 1:L_sensor2
    filename = sensor2{j}.name;
    found_i = false; % initialize a flag to indicate whether we've found i in the filename
    for i = 1:n
        if contains(filename, ['_Lt_t' num2str(i)])
            results_sensor2(i) = 1; % set the ith element of the results vector to 1
            found_i = true; % set the flag to true
            k=i;
            break; % break out of the inner loop once we've found i

        elseif contains(filename, ['_Rt_t' num2str(i)])
            results_sensor2(i+5) = 1; % set the ith element of the results vector to 1
            found_i = true; % set the flag to true
            k=i;
            break; % break out of the inner loop once we've found i

        elseif i>k && contains(filename, ['_Lt_t' num2str(i)])==0
            results_sensor2(i) = 0; % set the ith element of the results vector to 1

        elseif i>k && contains(filename, ['_Rt_t' num2str(i)])==0
            results_sensor2(i+5) = 0; % set the ith element of the results vector to 1
        end
    end
    if ~found_i
        disp(['The filename "', filename, '" does not contain any of the numbers 1 to ', num2str(n)])
    end
    if all(results_sensor2) % check if all elements of results are nonzero
        disp('All desired numbers have been found!')
        break; % break out of the outer loop if we've found all the numbers
    end
end


%% Making the sensor 1 and 2 matrix as 10 cell long.
sensor1_data_long=[sensor1; cell((L_drift-L_sensor1),1)];
k=0;
for i=1:10
    if results_sensor1(:,i)==1
        k=k+1;
    end
    if results_sensor1(:,i)==1 && (size(sensor1{k}.signals,2)>(critical_n_MU)); % discarding the trials when the detected no. of MU is less than critical value
        sensor1_data_long{i,:}=sensor1{k};
    else
        sensor1_data_long{i,:}=0;
        results_sensor1(:,i)=0;
    end
end

sensor2_data_long=[sensor2; cell((L_drift-L_sensor2),1)];
k=0;
for j=1:10
    if results_sensor2(:,j)==1
        k=k+1;
    end
    if results_sensor2(:,j)==1 && (size(sensor2{k}.signals,2)>(critical_n_MU)); % discarding the trials when the detected no. of MU is less than critical value
        sensor2_data_long{j,:}=sensor2{k};
        
    else
        sensor2_data_long{j}=0;
        results_sensor2(:,j)=0;
    end
end

%% emptying the force trial data wrt sensor1 and sensor2
null_sensor1=find(~results_sensor1);
null_sensor2=find(~results_sensor2);
null_sensor=unique([null_sensor1,null_sensor2]);
L_match=L_drift-size(null_sensor,2);
match_F_sensor=drift_f;
match_sensor1=sensor1_data_long;
match_sensor2=sensor2_data_long;
match_F_sensor(null_sensor)=[];
match_sensor1(null_sensor)=[];
match_sensor2(null_sensor)=[];