clc
clear
close all

%% Dynamical system
parameters.K=1;
parameters.M=1;
parameters.b=1;
parameters.Ts=0.1;
K=parameters.K;
M=parameters.M;
b=parameters.b;
Ts=parameters.Ts;
A=[1 Ts; -K*Ts/M 1-b*Ts/M];
B=[0; Ts/M];
clear K M b Ts
[parameters.n,parameters.m]=size(B);

%% LQR
parameters.Q=eye(parameters.n);
parameters.R=eye(parameters.M);
parameters.S=eye(parameters.n);

N=50*2;

x_0 = [8; -25];

oc1 = OC_util(parameters.Q, parameters.R, zeros(parameters.m, parameters.n), parameters.S, A, B);

[oc1.PP, oc1.FF] = OC_finite_horizon(oc1,N);

[x,u] = simulation(oc1, N, x_0, 0);

interval = 0:N-1;

figure; plot(interval,x(1,:),interval,x(2,:),interval,u); legend('state x1','state x2','input'); grid on; 
title('Plot States and Input'); xlabel('Time'); ylabel('Value');

