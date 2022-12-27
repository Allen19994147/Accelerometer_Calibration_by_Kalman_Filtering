%%% M271A Probability and Stochastic Process of Dyanmic Systems
%%% Allen Lee 705896702
%%% One realization
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
run_time = 10;%should be 30
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

%%% delta states IC
delta_states_accmeter = [x0;v0;bias_bar];
delta_states_accmeter_bar = zeros(3,1);
estimated_states = current_IMU_states;
M = [Mx0 0 0;0 Mv0 0;0 0 Mbias];
counter = 0;
err = zeros(2,1);
est_err_all = zeros(2,round(Num_Sample/(freq_acclerometer/freq_GPS))+1);
dy_err_all = zeros(2,Num_Sample+1);% A prior expected err
Num_Ensemble = 2000;
xk = zeros(2,1);
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

        %%% Expected values of a priori errors
        xk = current_true_states(1:2,1)...
            - (current_IMU_states(1:2,1) + delta_states_accmeter_bar(1:2,1));
        err = xk*xk';
        dy_err_all(1,k+1) = dy_err_all(1,k+1) + err(1,1);
        dy_err_all(2,k+1) = dy_err_all(2,k+1) + err(2,2);

        if(rem(k,(freq_acclerometer/freq_GPS))==0)% When measurements come...
            counter = counter + 1;
            Z_measurements = current_true_states(1:2,1)...
                +[normrnd(0,sqrt(Z_Variance(1,1)));normrnd(0,(sqrt(Z_Variance(2,2))))];
            delta_Z = Z_measurements - current_IMU_states(1:2,1);
            %%%%% Kalman algorithm part   %%%%%
            Kalman_Gain = M*H'/(H*M*H'+Z_Variance);
            P = M - (M*H'/(H*M*H'+Z_Variance))*H*M;
            M = P;
            delta_states_accmeter = delta_states_accmeter_bar...
                + Kalman_Gain*(delta_Z - H*delta_states_accmeter_bar);
            estimated_states(1:2,1) = current_IMU_states(1:2,1) + delta_states_accmeter(1:2,1);
            estimated_states(3,1) = delta_states_accmeter(3,1);
            %%% Ensemble check
            xk = current_true_states(1:2,1) - estimated_states(1:2,1);

            err = xk*estimated_states(1:2,1)';
            est_err_all(1,counter) = est_err_all(1,counter) + err(1,1);
            est_err_all(2,counter) = est_err_all(2,counter) + err(2,2);
            
        else
            delta_states_accmeter = delta_states_accmeter_bar;
            %%% Update estimated states %%%
            estimated_states(1:2,1) = current_IMU_states(1:2,1) + delta_states_accmeter(1:2,1);
            estimated_states(3,1) = delta_states_accmeter(3,1);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    est_err_all = est_err_all./Num_Ensemble;% Expected Values
    dy_err_all = dy_err_all./(Num_Sample+1);
end

%%

figure (2)
hold on
% plot(est_err_all(1,:),"red")
% plot(est_err_all(2,:),"blue")
% plot(dy_err_all(1,:),"red")
% plot(dy_err_all(2,:),"blue")
hold off
t1 = '$E[(x_k-H\bar{x_k})(x_j-H\bar{x_j})^{T}]$';
t2 = '$E[(x_k-\hat{x_k})\hat{x_k}^{T}]$';
title(t2,'interpreter','latex')
legend("Position Error","Velocity Error")
xlabel("Measurement Sample")
ylabel("Errors")

