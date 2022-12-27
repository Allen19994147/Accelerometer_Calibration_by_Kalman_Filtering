%%% M271A Probability and Stochastic Process of Dyanmic Systems
%%% Allen Lee 705896702
%%% More orthogonalities check
clc
close all
%%%%    True Model  %%%%%
s = rng; %random seed
x0_bar = 0;
Mx0 = 10^2;%100
v0_bar = 100;%100
Mv0 = 1;%1
bias_bar = 0;
Mbias = 0.01;

x0 = normrnd(x0_bar,sqrt(Mx0)); % Initial position
v0 = normrnd(v0_bar,sqrt(Mv0)); % Initial velocity
bias = normrnd(bias_bar,sqrt(Mbias));

% Variances of measurement(Z) noises
Z_Variance = [1 0;0 0.04^2];%1,0.04^2

w_bar = 0;
Mw = 0.0004;% 0.0004
%%% Parameters  %%%
A = 10;
omega = 0.1;% freq of true acc model 0.1
run_time = 30;%should be 30
freq_GPS = 5;% Hz, frequency of measurement
freq_acclerometer = 200;%Hz
delta_t = 1/freq_acclerometer;
Num_Sample = round(run_time*freq_acclerometer)+1;
time = (0:Num_Sample)*delta_t;

% True transistion Matrix
State_Transit = [1 delta_t 0.5*delta_t^2;0 1 delta_t;0 0 1];
% Accelerometer transition matrix
State_Transit_Dynamic = [1 delta_t -0.5*delta_t^2;0 1 -delta_t;...
    0 0 1];
Accelerometer_Noise_Matrix = -1.*[0.5*delta_t^2 delta_t 0]';
H = [1 0 0;0 1 0];% Measurement_Matrix
%%% Variables   %%%
current_true_states = [x0;v0;0];%p,v,a
current_IMU_states = [x0_bar;v0_bar;0];% Should be this one!

%%% Variables
Num_Ensemble = 1000;
delta_states_accmeter = [x0;v0;bias_bar];
delta_states_accmeter_bar = zeros(3,1);
estimated_states = current_IMU_states;
M = [Mx0 0 0;0 Mv0 0;0 0 Mbias];
counter = 0;
e_hat = zeros(2,1);
e_bar = zeros(2,1);
x_bar = zeros(2,1);
x_hat = zeros(2,1);
err_hat = zeros(3,1);% For ensemble variance P
err_hat_all = zeros(3*Num_Ensemble,round(Num_Sample/(freq_acclerometer/freq_GPS))+1);
P_all = zeros(3,round(Num_Sample/(freq_acclerometer/freq_GPS))+1); % P_all compares to err_hat_all
est_err_all = zeros(6,round(Num_Sample/(freq_acclerometer/freq_GPS))+1);
dy_err_all = zeros(4,Num_Sample+1);% A prior expected err
x_hat_all = zeros(3,round(Num_Sample/(freq_acclerometer/freq_GPS))+1);% for the last two orthogonal check
residual = zeros(2*Num_Ensemble,round(Num_Sample/(freq_acclerometer/freq_GPS))+1); % for residula term checking
for nn = 1 : Num_Ensemble   % # ensemble
    current_true_states = [x0;v0;0];%p,v,a
    current_IMU_states = [x0_bar;v0_bar;0];% Should be this one!
    delta_states_accmeter = [x0;v0;bias_bar];
    delta_states_accmeter_bar = zeros(3,1);
    estimated_states = current_IMU_states;
    M = [Mx0 0 0;0 Mv0 0;0 0 Mbias];
    err = zeros(2,1);
    counter = 0;

    for k = 0:Num_Sample
        %%% Calculate True model %%%
        current_true_states(3,1) = A*sin(k*delta_t*omega);% True acc
        current_true_states = State_Transit*current_true_states;
        %%% Calculate IMU Accelerometer  %%%
        w = normrnd(w_bar,sqrt(Mw));
        current_IMU_states(3,1) = A*sin(k*delta_t*omega)+bias+w;% IMU acc
        current_IMU_states = State_Transit*current_IMU_states;
        % A prior propogation of states(diff bt true/IMU)
        delta_states_accmeter_bar = State_Transit_Dynamic*delta_states_accmeter;

        M = State_Transit_Dynamic*M*State_Transit_Dynamic'...
            + Accelerometer_Noise_Matrix*Mw*Accelerometer_Noise_Matrix';

        %%% Expected values of a priori errors and...
        x_bar = delta_states_accmeter_bar(1:2,1);
        e_bar = current_true_states(1:2,1)-current_IMU_states(1:2,1) - x_bar;
        err = e_bar*e_bar';
        dy_err_all(1,k+1) = dy_err_all(1,k+1) + err(1,1);
        dy_err_all(2,k+1) = dy_err_all(2,k+1) + err(2,2);

        err = e_bar*x_bar';
        dy_err_all(3,k+1) = dy_err_all(3,k+1) + err(1,1);
        dy_err_all(4,k+1) = dy_err_all(4,k+1) + err(2,2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(rem(k,(freq_acclerometer/freq_GPS))==0)% When measurements come...
            counter = counter + 1;
            Z_measurements = current_true_states(1:2,1)...
                +[normrnd(0,sqrt(Z_Variance(1,1)));normrnd(0,(sqrt(Z_Variance(2,2))))];
            delta_Z = Z_measurements - current_IMU_states(1:2,1);
            %%% Residual check
            residual(2*nn-1:2*nn,counter) = delta_Z - x_bar;
            %%%%% Kalman algorithm part   %%%%%
            Kalman_Gain = M*H'/(H*M*H'+Z_Variance);
            P = M - (M*H'/(H*M*H'+Z_Variance))*H*M;
            M = P;
            P_all(:,3*counter-2:3*counter) = P;
            delta_states_accmeter = delta_states_accmeter_bar...
                + Kalman_Gain*(delta_Z - H*delta_states_accmeter_bar);
            x_hat_all(:,counter) = delta_states_accmeter;%% For last 2 ortho check
            estimated_states(1:2,1) = current_IMU_states(1:2,1) + delta_states_accmeter(1:2,1);
            estimated_states(3,1) = delta_states_accmeter(3,1);
            %%% Ensemble check  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            x_hat = delta_states_accmeter(1:2,1);
            e_hat = current_true_states(1:2,1)-current_IMU_states(1:2,1) - x_hat;
            err_hat = [e_hat;bias-estimated_states(3,1)];
            err_hat_all(3*nn-2:3*nn,counter) = err_hat;
            err = e_hat*x_hat';
            est_err_all(1,counter) = est_err_all(1,counter) + err(1,1);
            est_err_all(2,counter) = est_err_all(2,counter) + err(2,2);
            err = e_hat*delta_Z';
            est_err_all(3,counter) = est_err_all(3,counter) + err(1,1);
            est_err_all(4,counter) = est_err_all(4,counter) + err(2,2);
            err = e_hat*x_bar';
            est_err_all(5,counter) = est_err_all(5,counter) + err(1,1);
            est_err_all(6,counter) = est_err_all(6,counter) + err(2,2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            delta_states_accmeter = delta_states_accmeter_bar;
            %%% Update estimated states %%%
            estimated_states(1:2,1) = current_IMU_states(1:2,1) + delta_states_accmeter(1:2,1);
            estimated_states(3,1) = delta_states_accmeter(3,1);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

end % End resemble
est_err_all = est_err_all./Num_Ensemble;% Expected Values
dy_err_all = dy_err_all./Num_Ensemble;
e_ave = zeros(3,length(est_err_all));
for i=1:counter % for each measurement
    for j = 1:Num_Ensemble% averaging over it
        e_ave(:,i) = e_ave(:,i) + err_hat_all(3*j-2:3*j,i);
    end
    e_ave(:,i) = e_ave(:,i)./Num_Ensemble;
end
P_esmb = zeros(3,3*counter);
x_telda = zeros(3,1);
for i=1:length(e_ave)% for each measurement
    for j = 1:length(err_hat_all)/3 % 1 to 3*Num_Ensemble
        x_telda = err_hat_all(3*j-2:3*j,i) - e_ave(:,i);
        P_esmb(:,3*i-2:3*i) = P_esmb(:,3*i-2:3*i) + x_telda*x_telda';
    end
end
P_esmb = P_esmb./(Num_Ensemble-1);
%%

figure (2)
hold on
%  plot(est_err_all(1,:),"red")   %   e_hat*x_hat'
% plot(est_err_all(3,:),"red")  %   e_hat*delta_Z'
% plot(est_err_all(5,:),"red")  %   e_hat*x_bar'

plot(est_err_all(2,:),"blue")
% plot(est_err_all(4,:),"blue")
% plot(est_err_all(6,:),"blue")

% plot(dy_err_all(1,:),"red")   % useless, don't use % e_bar*e_bar'
%plot(dy_err_all(3,:),"red")   % e_bar*x_bar'

% plot(dy_err_all(2,200:end),"blue")
% plot(dy_err_all(4,200:end),"blue")

hold off
title("$E[\hat e_k \delta Z_k^T],Velocity$",'Interpreter','latex')
ylabel("Expected Value",'Interpreter','latex')
xlabel("Measurement time step",'Interpreter','latex')

%%
%%% Compare P and P_esmb
P_check = zeros(3,counter);
P_null = P_all - P_esmb;
for i=1:length(P_all)/3 % Extract the diagonal terms
    for j = 1:3
        P_check(j,i) = P_null(j,3*i-3+j);
    end
end
figure(3)
hold on
% plot(P_check(1,:),"red")
% plot(P_check(2,:),"blue")
plot(P_check(3,:),"green")
hold off
title("Error of bias variance")
ylabel("Bias error $(m/s^2)^2$",'Interpreter','latex')
xlabel("Measurement time step",'Interpreter','latex')

%%
%%% Compare the last two orthogonality check
lasttwocheck_long = zeros(3,counter*3);
for i=1:length(e_ave)% for each measurement
    for j = 1:length(err_hat_all)/3 % 1 to 3*Num_Ensemble
        x_telda = err_hat_all(3*j-2:3*j,i) - e_ave(:,i);
        lasttwocheck_long(:,3*i-2:3*i) = lasttwocheck_long(:,3*i-2:3*i) + x_telda*x_hat_all(:,i)';
    end
end
lasttwocheck_long = lasttwocheck_long./(Num_Ensemble-1);
lasttwocheck = zeros(3,counter);
for i=1:length(P_all)/3 % Extract the diagonal terms
    for j = 1:3
        lasttwocheck(j,i) = lasttwocheck_long(j,3*i-3+j);
    end
end
figure(4)
hold on
% plot(lasttwocheck(1,:),"red")
% plot(lasttwocheck(2,:),"blue")
plot(lasttwocheck(3,:),"green")
hold off
title("Bias Orthogonality")
ylabel("Error",'Interpreter','latex')
xlabel("Measurement time step",'Interpreter','latex')
%%
%%% Perform the residual check
% We use k-1 and k to verify
Residual_variance = zeros(2,2*counter-2);
for i=2:counter
    for j = 1:Num_Ensemble
        Residual_variance(:,2*i-3:2*i-2) = residual(2*j-1:2*j,i-1)*residual(2*j-1:2*j,i)';
    end
end
Residual_variance = Residual_variance./Num_Ensemble;
Resi_lastcheck = zeros(2,counter-1);
for i=1:length(Residual_variance)/2 % Extract the diagonal terms
    for j = 1:2
        Resi_lastcheck(j,i) = Residual_variance(j,2*i-2+j);
    end
end

figure(5)
hold on
plot(Resi_lastcheck(1,:),"red")
% plot(Resi_lastcheck(2,:),"blue")
hold off
title("Position of the residual")
ylabel("Residual error",'Interpreter','latex')
xlabel("Measurement time step",'Interpreter','latex')