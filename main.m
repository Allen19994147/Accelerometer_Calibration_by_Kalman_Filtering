%%% M271A Probability and Stochastic Process of Dyanmic Systems
%%% Allen Lee 705896702
%%% Toy
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

x0 = normrnd(x0_bar,Mx0); % Initial position
v0 = normrnd(v0_bar,Mv0); % Initial velocity
bias = normrnd(bias_bar,Mbias);

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
% current_IMU_states = current_true_states;% Change it!

%%% delta states IC
delta_states_accmeter = [x0;v0;bias_bar];
%delta_states_accmeter = [0;0;0];
delta_states_accmeter_bar = zeros(3,1);
estimated_states = current_IMU_states;
M = [Mx0 0 0;0 Mv0 0;0 0 Mbias];
counter = 0;
error_all = zeros(3,Num_Sample+1);
err = 0;
sigma_bound = 0;
simga_bound_all = zeros(3,Num_Sample+1);
filter_err_var = zeros(3,Num_Sample+1);

figure(1)
hold on
for k = 0:Num_Sample
    %%% Calculate True model %%%
    current_true_states(3,1) = A*sin(k*delta_t*omega);% True acc
    current_true_states = State_Transit*current_true_states;
    plot(current_true_states(1),current_true_states(2),".",Color="blue")

    %%% Calculate IMU Accelerometer  %%%
    w = normrnd(w_bar,Mw);
    current_IMU_states(3,1) = A*sin(k*delta_t*omega)+bias+w;% IMU acc
    current_IMU_states = State_Transit*current_IMU_states;

    % A prior propogation of states(diff bt true/IMU)
    delta_states_accmeter_bar = State_Transit_Dynamic*delta_states_accmeter;

    M = State_Transit_Dynamic*M*State_Transit_Dynamic'...
        + Accelerometer_Noise_Matrix*Mw*Accelerometer_Noise_Matrix';

    if(rem(k,(freq_acclerometer/freq_GPS))==0)% When measurements come...
        Z_measurements = current_true_states(1:2,1)...
            +[normrnd(0,Z_Variance(1,1));normrnd(0,Z_Variance(2,2))];
        delta_Z = Z_measurements - current_IMU_states(1:2,1);
        plot(Z_measurements(1),Z_measurements(2),".",Color="green")
%%%%% Kalman algorithm part   %%%%%
        Kalman_Gain = M*H'/(H*M*H'+Z_Variance);
        P = M - (M*H'/(H*M*H'+Z_Variance))*H*M;
        M = P;
        delta_states_accmeter = delta_states_accmeter_bar...
            + Kalman_Gain*(delta_Z - H*delta_states_accmeter_bar);

        estimated_states(1:2,1) = current_IMU_states(1:2,1) + delta_states_accmeter(1:2,1);
        estimated_states(3,1) = delta_states_accmeter(3,1);

        plot(estimated_states(1),estimated_states(2),".",Color="red",MarkerSize=10)

    else
        delta_states_accmeter = delta_states_accmeter_bar;
        %%% Update estimated states %%%
        estimated_states(1:2,1) = current_IMU_states(1:2,1) + delta_states_accmeter(1:2,1);
        estimated_states(3,1) = delta_states_accmeter(3,1);
        %plot(estimated_states(1),estimated_states(2),".",Color="yellow")
    end

w

    %%%%%%%%%%%     Error calculation   %%%%%%%%%%%%%%%%%
    err = (current_true_states-estimated_states);
    error_all(3,k+1) = bias - delta_states_accmeter(3,1);
    for i = 1:2
        error_all(i,k+1) = (err(i));
    end

    for i = 1:3
        simga_bound_all(i,k+1) = sqrt(M(i,i));
        filter_err_var(i,k+1) = M(i,i);
    end

end
hold off
legend("True states","Measurement","Estimated States")
title("GPS at 100Hz")
xlabel("Position (m)")
ylabel("Velocity (m/s)")

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
title("Error of bias")
xlabel("Time step $\bigtriangleup t$",'Interpreter','latex')
ylabel("Bias Error $(m/s^2)$",'Interpreter','latex')

