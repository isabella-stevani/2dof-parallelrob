% Author: Isabella Stevani, Masters Student
% Master's Dissertation Project: Robust Control of a 2-DOF Parallel
% Mechanism Combining Feedback Linearization and H-Infinity Design
% Code: H-Infinity Design
% Polytechnic School of The University of Sao Paulo, Dept. of 
% Telecommunications and Control (PTC)
% E-mail address: isabella.stevani@usp.br
% Creation: Dec 2019; Last revision: 10-Dec-2019

close all; clear; clc;

%------------- BEGIN OF CODE --------------

%% Linear estimation for uncertain model

% Model parameters

%%% Nominal parameters (SI units)
g = 9.8;
m1 = 0.1;
m2 = 0.1;
l0 = 0.05;
l1 = 0.12;
l2 = 0.16;
lg1 = 0.06;
lg2 = 0.068;
Jz1 = 2e-04;
Jz2 = 7e-04;

param.g = g;
param.m1 = m1;
param.m2 = m2;
param.l0 = l0;
param.l1 = l1;
param.l2 = l2;
param.lg1 = lg1;
param.lg2 = lg2;
param.Jz1 = Jz1;
param.Jz2 = Jz2;

lambda = 40; %FL parameter [rad/s]

% Simulation parameters

T = 0.001; %sample time [s]
tsim = 1; %simulation time [s]
t = 0:T:tsim; %simulation time vector
sat = 2; %actuators' saturation

%%% Initial conditions
r_xy = [0.07; 0.17]; %end-effector initial condition
rad72 = 72*pi/180;
q12_0 = [rad72; rad72; rad72; rad72]; %joint angles initial conditions
inputfunc = @ConstantInput; %constant reference signal for initial conditions

%%% Adjusting initial conditions through kinematics
q = ParallelRobKinematics(r_xy,q12_0,param,inputfunc,lambda);
dq = zeros(6,1);
x0 = [q;dq]; %initial state vector

% Input data

%%% Parameters
cord = 'theta'; %reference considering end-effector coordinates
inputfunc = @RoundInput; %round reference signal
% inputfunc = @ConstantInput; %constant reference signal

%%% Signal computation
ref = ParallelRobDynRef(x0,t,param,inputfunc,lambda);

% Output data

loop = 'OL'; %open-loop simulation
FF.on = true; %with feedforward

%%% Dynamics
[q,~,~] = ParallelRobDynamics(x0,t,param,param,lambda, ...
    ref,cord,T,sat,FF,loop);

% Transfer function

%%% Data
if strcmp(cord,'xy')
    Qa = [1 0;
        0 1;
        0 0;
        0 0;
        0 0;
        0 0];
elseif strcmp(cord,'theta')
    Qa = [0 0;
        0 0;
        1 0;
        0 0;
        0 1;
        0 0];
end
in = ref;
out = q;
%% Nominal model

% Model parameters
lambda = 40; %FL parameter [rad/s]
s = tf('s');

% Continuous linearized model
G = 1/(s+lambda)^2;

figure;
sigma(G);
%------------- END OF CODE --------------