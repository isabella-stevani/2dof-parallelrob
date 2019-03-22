% Author: Isabella Stevani, Masters Student
% Master's Dissertation Project: Robust Control of a 2-DOF Parallel
% Mechanism Combining Feedback Linearization and H-Infinity Design
% Code: Discrete Control Design and Simulationclos
% Polytechnic School of The University of Sao Paulo, Dept. of 
% Telecommunications and Control (PTC)
% E-mail address: isabella.stevani@usp.br
% Creation: Aug 2018; Last revision: 18-Mar-2019

close all; clear; clc;

%------------- BEGIN OF CODE --------------

%% Folder paths

% Determine where your m-file's folder is
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path
addpath(genpath(folder));

%% Parallel model parameters

T = 0.001; %sample time [s]

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

%%% Nominal linearized model
lambda = 40; %FL parameter [rad/s]
z = tf('z');

% Discrete model
z1 = 0.82; z2 = 0.8; z3 = 0.22;
Gd = (T^2/8)*((1-z1)*(1-z2)*(1-z3)*z*(1+z)^3)/((-1+z)^2*(-z1+z)*(-z2+z)*(-z3+z));
Gd.Ts = T;

%% Simulation parameters

tsim = 1; %simulation time [s]
t = 0:T:tsim; %simulation time vector
sat = 2; %actuators' saturation

%%% Initial conditions
r_xy = [0.07; 0.17]; %end-effector initial condition
rad72 = 72*pi/180;
q12_0 = [rad72; rad72; rad72; rad72]; %joint angles initial conditions
inputfunc = @ConstantInput; %constant reference signal for initial conditions

% Adjusting initial conditions through kinematics
q = ParallelRobKinematics(r_xy,q12_0,param,inputfunc,lambda);
dq = zeros(6,1);
x0 = [q;dq]; %initial state vector

%%% Reference signal
% Parameters
cord = 'theta'; %reference considering end-effector coordinates
% cord = 'xy'; %reference considering active joints coordinates
inputfunc = @RoundInput; %round reference signal
% inputfunc = @ConstantInput; %constant reference signal

% Signal computation
ref = ParallelRobDynRef(x0,t,param,inputfunc,lambda,cord);
% save('ref_theta.mat','ref');
% load('ref_theta');

%% Nominal simulation

%%% Dynamics
[q,dq,u] = ParallelRobDynamics_d(x0,t,param,param,lambda,ref,cord,T,sat);

% States
figure; leg = {'FL'};
ylab = {'x [m]','y [m]','t11 [rad]','t12 [rad]','t21 [rad]','t22 [rad]'};
xlab = 'Time [s]';
tlab = 'States';
for i = 1:6
    subplot(3,2,i)
    plot(t,q(i,:),'-b');
    grid on;
%     lgd = legend(leg);
    xl = xlabel(xlab);
    yl = ylabel(ylab{i});
    set(xl,'FontSize',16);
    set(yl,'FontSize',16);
%     set(lgd,'FontSize',12);
end
supertitle(tlab,18);
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 2;
end
set(gcf,'PaperPositionMode','auto');
saveas(gcf,['images/cntr_nomsim_states','.eps'], 'epsc');
saveas(gcf,'images/cntr_nomsim_states.png');
% 
% Control signal
figure; leg = {'FL'};
xlab = {'','Time [s]'};
ylab = {'u1 [N.m]','u2 [N.m]'};
tlab = 'Control Signals';
for i = 1:2
    subplot(2,1,i)
    plot(t,u(i,:),'-b');
    grid on;
%     lgd = legend(leg);
    xl = xlabel(xlab{i});
    yl = ylabel(ylab{i});
%     set(lgd,'FontSize',12);
    set(xl,'FontSize',16);
    set(yl,'FontSize',16);
end
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 2;
end
supertitle(tlab,18);
set(gcf,'PaperPositionMode','auto');
saveas(gcf,['images/cntr_nomsim_cntrsig','.eps'], 'epsc');
saveas(gcf,'images/cntr_nomsim_cntrsig.png');

% Actuated states errors
r = ref(1:2,:);
e = r-[q(3,:);q(5,:)];
figure; leg = {'FL'};
xlab = {'','Time [s]'};
ylab = {'e1 [rad]','e2 [rad]'};
tlab = 'Actuated Angle Errors';
for i = 1:2
    subplot(2,1,i)
    plot(t,e(i,:),'-b');
    grid on;
%     lgd = legend(leg);
    xl = xlabel(xlab{i});
    yl = ylabel(ylab{i});
%     set(lgd,'FontSize',12);
    set(xl,'FontSize',16);
    set(yl,'FontSize',16);
end
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 2;
end
supertitle(tlab,18);
set(gcf,'PaperPositionMode','auto');
saveas(gcf,['images/cntr_nomsim_error','.eps'], 'epsc');
saveas(gcf,'images/cntr_nomsim_error.png');

%------------- END OF CODE --------------