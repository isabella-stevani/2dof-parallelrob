% Author: Isabella Stevani, Masters Student
% Master's Dissertation Project: Robust Control of a 2-DOF Parallel
% Mechanism Combining Feedback Linearization and H-Infinity Design
% Code: H-Infinity Design
% Polytechnic School of The University of Sao Paulo, Dept. of 
% Telecommunications and Control (PTC)
% E-mail address: isabella.stevani@usp.br
% Creation: Dec 2019; Last revision: 07-Jan-2020

close all; clear; clc;

%------------- BEGIN OF CODE --------------

%% Parameters

set_env; %script to set work environment

%% Nominal simulation

tsim = 1; %simulation time [s]
t = 0:T:tsim; %simulation time vector

% Reference signal

ref = ParallelRobDynRef(x0,t,param,inputfunc,lambda);

% Dynamics

x0 = x0+x0_tol; %initial conditions tolerance
sim_type = 'design'; %control design simulation
FF.on = false; %without feed-forward

[q,~,~] = ParallelRobDynamics(x0,t,param,param,lambda, ...
    ref,cord,T,sat,FF,sim_type);

% Transfer function estimation

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

in = ref(1:6,:);
out = q;

in_a = zeros(2,length(t));
out_a = zeros(2,length(t));

for i = 1:length(t)
    in_a(:,i) = Qa'*in(:,i);
    out_a(:,i) = Qa'*out(:,i);
end
in_a = in_a';
out_a = out_a';
x0_a = Qa'*x0(1:6);
data = iddata(out_a,in_a,T);
data1 = iddata(out_a(:,1),in_a(:,1),T);
data2 = iddata(out_a(:,2),in_a(:,2),T);

%%% Estimation
% tfest_opt = tfestOptions;
% tfest_opt.InitializeMethod = 'all';
sys = tfest(data,2,0);
sys1 = tfest(data1,2,0);
sys2 = tfest(data2,2,0);

%% Nominal model

% Model parameters
s = tf('s');

% Continuous linearized model
G_OL = ((lambda^2)/(s+lambda)^2)*eye(2);
% G_CL = feedback(G_OL,eye(2));
G_CL = G_OL;
G_CL_ss = ss(G_CL);
[y,~,~] = lsim(G_CL_ss,in_a,t,[x0_a(1);x0(9);x0_a(2);x0(11)]);

G_CL_est = [sys1 0;0 sys2];
G_CL_est_ss = ss(G_CL_est);
[y_est,~,~] = lsim(G_CL_est_ss,in_a,t,[x0_a(1);x0(9);x0_a(2);x0(11)]);

% figure;
% sigma(G,sys(1,1),sys(2,2));
% legend on;
% figure;
% sigma(G,sys1);
% legend on;

figure; hold on;
plot(t,y(:,1),'-b');
plot(t,y_est(:,1),'-k');
plot(t,out_a(:,1),'-r');
plot(t,in_a(:,1),'--k');
legend('TF','TF est','model','ref');
hold off;

figure; hold on;
plot(t,y(:,2),'-b');%%
plot(t,y_est(:,2),'-k');
plot(t,out_a(:,2),'-r');
plot(t,in_a(:,2),'--k');
legend('TF','TF est','model','ref');
hold off;

% figure; hold on;
% plot(ref(1,:),ref(2,:),'--k');
% plot(out(1,:),out(2,:),'-r');
% hold off;
% 
% e1 = ref(1,:)-out(1,:);
% e2 = ref(2,:)-out(2,:);
% figure;
% subplot(2,1,1);
% plot(t,e1);
% subplot(2,1,2);
% plot(t,e2);

%% Uncertain model

tolI = 0.2; % 20% accurate (inertia parameters)
toll = 0.01; % 1% accurate (bar lengths)
p = 4;

uncparam.g = g;
uncparam.m1 = round(m1*(1+tolI),p);
uncparam.m2 = round(m2*(1+tolI),p);
uncparam.l0 = round(l0*(1+toll),p);
uncparam.l1 = round(l1*(1+toll),p);
uncparam.l2 = round(l2*(1+toll),p);
uncparam.lg1 = round(lg1*(1+tolI),p);
uncparam.lg2 = round(lg2*(1+tolI),p);
uncparam.Jz1 = round(Jz1*(1+tolI),p);
uncparam.Jz2 = round(Jz2*(1+tolI),p);

[q_unc,~,~] = ParallelRobDynamics(x0,t,param,uncparam,lambda, ...
    ref,cord,T,sat,FF,sim_type);

in = ref(1:6,:);
out = q_unc;

in_a = zeros(2,length(t));
out_a = zeros(2,length(t));

for i = 1:length(t)
    in_a(:,i) = Qa'*in(:,i);
    out_a(:,i) = Qa'*out(:,i);
end
in_a = in_a';
out_a = out_a';
x0_a = Qa'*x0(1:6);
data = iddata(out_a,in_a,T);
data1 = iddata(out_a(:,1),in_a(:,1),T);
data2 = iddata(out_a(:,2),in_a(:,2),T);

sys = tfest(data,2,0);
sys1 = tfest(data1,2,0);
sys2 = tfest(data2,2,0);

G_CL_est = [sys1 0;0 sys2];
G_CL_est_ss = ss(G_CL_est);
[y_est,~,~] = lsim(G_CL_est_ss,in_a,t,[x0_a(1);x0(9);x0_a(2);x0(11)]);

% figure;
% sigma(G,sys(1,1),sys(2,2));
% legend on;
% figure;
% sigma(G,sys1);
% legend on;

figure; hold on;
plot(t,y(:,1),'-b');
plot(t,y_est(:,1),'-k');
plot(t,out_a(:,1),'-r');
plot(t,in_a(:,1),'--k');
legend('TF','TF est','model','ref');
hold off;

figure; hold on;
plot(t,y(:,2),'-b');%%
plot(t,y_est(:,2),'-k');
plot(t,out_a(:,2),'-r');
plot(t,in_a(:,2),'--k');
legend('TF','TF est','model','ref');
hold off;

%------------- END OF CODE --------------