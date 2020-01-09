% Author: Isabella Stevani, Masters Student
% Master's Dissertation Project: Robust Control of a 2-DOF Parallel
% Mechanism Combining Feedback Linearization and H-Infinity Design
% Code: Set work environment
% Polytechnic School of The University of Sao Paulo, Dept. of 
% Telecommunications and Control (PTC)
% E-mail address: isabella.stevani@usp.br
% Creation: Jan 2020; Last revision: 07-Jan-2020

% close all; clear; clc;

%------------- BEGIN OF CODE --------------

%% Folder paths

% Determine where your m-file's folder is
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path
addpath(genpath(folder));

%% Nominal parameters (SI units)

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

% Feedback linearization
FL.p1 = -40;
FL.p2 = -40;
FL.K1 = -(FL.p1+FL.p2);
FL.K2 = FL.p1*FL.p2;
FL.lambdabar = 10*max(abs(FL.p1),abs(FL.p2)); %quick convergence for coupling equations


%% Simulation parameters

T = 0.001; %sample time [s]
sat = 2; %actuators' saturation

% Initial conditions

%%% Parameters
r_xy = [0.07; 0.17]; %end-effector initial condition
rad72 = 72*pi/180;
q12_0 = [rad72; rad72; rad72; rad72]; %joint angles initial conditions
inputfunc = @ConstantInput; %constant reference signal for initial conditions

%%% Adjusting initial conditions through kinematics
q = ParallelRobKinematics(r_xy,q12_0,param,inputfunc,FL);
dq = zeros(6,1);
x0 = [q;dq]; %initial state vector
x0_tol = [0.1;0.1;0.2;0.2;0.2;0.2;0;0;0;0;0;0];

% Input data

%%% Parameters
cord = 'theta'; %reference considering end-effector coordinates
inputfunc = @RoundInput; %round reference signal

%------------- END OF CODE --------------