%%% M271A Probability and Stochastic Process of Dyanmic Systems
%%% Allen Lee 705896702
%%% 2 measurement resources
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The task is to demonstrate that batch update and sequential update are
% the same for Kalman Filter when multiple measurements
% sources are available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
bias = normrnd(bias_bar,sqrt(Mbias)); % const bias

% Variances of measurement(Z) noises
Z1_Variance = [1 0;0 0.04^2];% Measurement 1
Z2_Variance = [2 0;0 0.02^2];% Measurement 2
batch_Z_Variance = [Z1_Variance zeros(2,2);zeros(2,2) Z2_Variance];
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
H = [1 0 0;0 1 0];% Measurement_Matrix, sequential
HH = [H;H]; % Measurement_Matrix, sequential
%%% Variables   %%%
current_true_states = [x0;v0;0];%p,v,a
current_IMU_states = [x0_bar;v0_bar;0];% Should be this one!

%%% Initialization
delta_states_accmeter_bh = [x0;v0;bias_bar];
delta_states_accmeter_sq = [x0;v0;bias_bar];
delta_states_accmeter_bar_bh = zeros(3,1);
delta_states_accmeter_bar_sq = zeros(3,1);
estimated_states_bh = current_IMU_states;% sequential estimation
estimated_states_sq = current_IMU_states;% batch estimation
M_sq = [Mx0 0 0;0 Mv0 0;0 0 Mbias];
M_bh = [Mx0 0 0;0 Mv0 0;0 0 Mbias];
error_all = zeros(3,Num_Sample+1);
err = 0;
sigma_bound = 0;
simga_bound_all = zeros(3,Num_Sample+1);
filter_err_var = zeros(3,Num_Sample+1);
P = zeros(2,2);
Z1_measurements = [0;0];
Z2_measurements = [0;0];
batch_Z = [Z1_measurements;Z2_measurements];
figure(1)
hold on
for k = 0:Num_Sample
    %%% Calculate True model %%%
    current_true_states(3,1) = A*sin(k*delta_t*omega);% True acc
    current_true_states = State_Transit*current_true_states;
    plot(current_true_states(1),current_true_states(2),".",Color="blue")

    %%% Calculate IMU Accelerometer  %%%
    w = normrnd(w_bar,sqrt(Mw));
    current_IMU_states(3,1) = A*sin(k*delta_t*omega)+bias+w;% IMU acc
    current_IMU_states = State_Transit*current_IMU_states;

    % A prior propogation of states(diff bt true/IMU)
    delta_states_accmeter_bar_sq = State_Transit_Dynamic*delta_states_accmeter_sq;
    delta_states_accmeter_bar_bh = State_Transit_Dynamic*delta_states_accmeter_bh;
    M_sq = State_Transit_Dynamic*M_sq*State_Transit_Dynamic'...
        + Accelerometer_Noise_Matrix*Mw*Accelerometer_Noise_Matrix';
    M_bh = State_Transit_Dynamic*M_bh*State_Transit_Dynamic'...
        + Accelerometer_Noise_Matrix*Mw*Accelerometer_Noise_Matrix';

    if(rem(k,(freq_acclerometer/freq_GPS))==0)% When measurements come in...
        Z1_measurements = current_true_states(1:2,1)...
            +[normrnd(0,sqrt(Z1_Variance(1,1)));normrnd(0,(sqrt(Z1_Variance(2,2))))];
        Z2_measurements = current_true_states(1:2,1)...
            +[normrnd(0,sqrt(Z2_Variance(1,1)));normrnd(0,(sqrt(Z2_Variance(2,2))))];
        delta_Z1 = Z1_measurements - current_IMU_states(1:2,1);
        delta_Z2 = Z2_measurements - current_IMU_states(1:2,1);

        %%%%% Kalman algorithm part,Sequential update   %%%%%
        Kalman_Gain = M_sq*H'/(H*M_sq*H'+Z1_Variance);
        P = M_sq - (M_sq*H'/(H*M_sq*H'+Z1_Variance))*H*M_sq;
        M_sq = P;
        delta_states_accmeter_sq = delta_states_accmeter_bar_sq...
            + Kalman_Gain*(delta_Z1 - H*delta_states_accmeter_bar_sq);

        Kalman_Gain = M_sq*H'/(H*M_sq*H'+Z2_Variance);
        P = M_sq - (M_sq*H'/(H*M_sq*H'+Z2_Variance))*H*M_sq;
        M_sq = P;
        delta_states_accmeter_sq = delta_states_accmeter_sq...
            + Kalman_Gain*(delta_Z2 - H*delta_states_accmeter_sq);

        estimated_states_sq(1:2,1) = current_IMU_states(1:2,1) + delta_states_accmeter_sq(1:2,1);
        estimated_states_sq(3,1) = delta_states_accmeter_sq(3,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Kalman algorithm part,batch update (all at once)  %%%%%
        batch_Z = [delta_Z1;delta_Z2];
        Kalman_Gain = M_bh*HH'/(HH*M_bh*HH'+batch_Z_Variance);
        P = M_bh - (M_bh*HH'/(HH*M_bh*HH'+batch_Z_Variance))*HH*M_bh;
        M_bh = P;
        delta_states_accmeter_bh = delta_states_accmeter_bar_bh...
            + Kalman_Gain*(batch_Z - HH*delta_states_accmeter_bar_bh);


        estimated_states_bh(1:2,1) = current_IMU_states(1:2,1) + delta_states_accmeter_bh(1:2,1);
        estimated_states_bh(3,1) = delta_states_accmeter_bh(3,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Check if two update are the same %%%%%%%%%
        if(sum(abs(delta_states_accmeter_sq-delta_states_accmeter_bh))>0.01)
            delta_states_accmeter_sq - delta_states_accmeter_bh
            error("Sequential update and batch update not consistenet!\n")
        end
        plot(estimated_states_sq(1),estimated_states_sq(2),".",Color="red",MarkerSize=5)
        plot(estimated_states_bh(1),estimated_states_bh(2),"o",Color="green",MarkerSize=5)
    else
        delta_states_accmeter_sq = delta_states_accmeter_bar_sq;
        delta_states_accmeter_bh = delta_states_accmeter_bar_bh;
        %%% Update estimated states, sequential %%%
        estimated_states_sq(1:2,1) = current_IMU_states(1:2,1) + delta_states_accmeter_sq(1:2,1);
        estimated_states_sq(3,1) = delta_states_accmeter_sq(3,1);
        %%% Update estimated states, batch %%%
        estimated_states_bh(1:2,1) = current_IMU_states(1:2,1) + delta_states_accmeter_bh(1:2,1);
        estimated_states_bh(3,1) = delta_states_accmeter_bh(3,1);
    end



    %%%%%%%%%%%     Error calculation (batch)   %%%%%%%%%%%%%%%%%
    err = (current_true_states-estimated_states_bh);
    error_all(3,k+1) = bias - delta_states_accmeter_bh(3,1);
    for i = 1:2
        error_all(i,k+1) = (err(i));
    end

    for i = 1:3
        simga_bound_all(i,k+1) = sqrt(M_bh(i,i));
        filter_err_var(i,k+1) = M_bh(i,i);
    end

end
hold off
legend("True states","Seq Estimation","Batch Estimation")
title("True Model, and Estimation")
xlabel("Position (m)")
ylabel("Velocity (m/s)")
hold off
%%

figure (2)
hold on

plot(error_all(1,:),"red")
plot(simga_bound_all(1,:),"yellow")
plot(filter_err_var(1,:),"black")
plot(-1.*filter_err_var(1,:),"black")
plot(-1.*simga_bound_all(1,:),"yellow")

%
% plot(error_all(2,:),"blue")
% plot(-1.*simga_bound_all(2,:),"yellow")
% plot(filter_err_var(2,:),"black")
% plot(-1.*filter_err_var(2,:),"black")
% plot(simga_bound_all(2,:),"yellow")

%
% plot(error_all(3,:),"green")
% plot(simga_bound_all(3,:),"yellow")
% plot(filter_err_var(3,:),"black")
% plot(-1.*filter_err_var(3,:),"black")
% plot(-1.*simga_bound_all(3,:),"yellow")

hold off
legend("Error","1 sigma bound","Filter error variance")
title("Error of position")
xlabel("Time step $\bigtriangleup t$",'Interpreter','latex')
ylabel("Position Error $(m/s^2)$",'Interpreter','latex')

