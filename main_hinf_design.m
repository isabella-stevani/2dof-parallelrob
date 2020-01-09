% Author: Isabella Stevani, Masters Student
% Master's Dissertation Project: Robust Control of a 2-DOF Parallel
% Mechanism Combining Feedback Linearization and H-Infinity Design
% Code: H-Infinity Design
% Polytechnic School of The University of Sao Paulo, Dept. of 
% Telecommunications and Control (PTC)
% E-mail address: isabella.stevani@usp.br
% Creation: Dec 2019; Last revision: 09-Jan-2020

close all; clear; clc;

%------------- BEGIN OF CODE --------------

%% Parameters

set_env; %script to set work environment

%% Nominal case

% Simulation

tsim = 1; %simulation time [s]
t = 0:T:tsim; %simulation time vector

%%% Reference signal
ref = ParallelRobDynRef(x0,t,param,inputfunc,FL);

%%% Dynamics
x0 = x0+x0_tol; %initial conditions tolerance
sim_type = 'design'; %control design simulation
FF.on = false; %without feed-forward

[q,~,~] = ParallelRobDynamics(x0,t,param,param,FL, ...
    ref,cord,T,sat,FF,sim_type);

% Transfer function estimation

%%% Raw data
in = ref(1:6,:);
out = q;

%%% Actuated data
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

in_a = zeros(2,length(t));
out_a = zeros(2,length(t));

for i = 1:length(t)
    in_a(:,i) = Qa'*in(:,i);
    out_a(:,i) = Qa'*out(:,i);
end
in_a = in_a';
out_a = out_a';

x0_a = Qa'*x0(1:6); %actuated initial conditions

e_a = in_a-out_a;

% data = iddata(out_a,in_a,T);
data1 = iddata(out_a(:,1),in_a(:,1),T); %first actuated variable data
data2 = iddata(out_a(:,2),in_a(:,2),T); %second actuated variable data

%%% Estimation
% tfest_opt = tfestOptions;
% tfest_opt.InitializeMethod = 'all';
% sys = tfest(data,2,0);
sys1_nom = tfest(data1,2,0); %estimated TF for the first actuated variable
sys2_nom = tfest(data2,2,0); %estimated TF for the second actuated variable

%%% Estimated TF response
G_est_nom = [sys1_nom 0;0 sys2_nom];
G_est_nom_ss = ss(G_est_nom);
[y_est_nom,~,~] = lsim(G_est_nom_ss,in_a,t,[x0_a(1);x0(9);x0_a(2); ...
    x0(11)]);

% Nominal model

%%% Continuous linearized model
s = tf('s');
G = FL.K2/(s^2+FL.K1*s+FL.K2)*eye(2);
G_ss = ss(G);
[y,~,~] = lsim(G_ss,in_a,t,[x0_a(1);x0(9);x0_a(2);x0(11)]);

% Plots

% figure;
% sigma(G,sys(1,1),sys(2,2));
% legend on;
% figure;
% sigma(G,sys1);
% legend on;

figure; hold on;
% plot(t,y(:,1),'-b');
plot(t,y_est_nom(:,1),'-k');
plot(t,out_a(:,1),'-r');
plot(t,in_a(:,1),'--k');
% legend('TF','TF est','model','ref');
legend('TF est','model','ref');
hold off;

figure; hold on;
% plot(t,y(:,2),'-b');
plot(t,y_est_nom(:,2),'-k');
plot(t,out_a(:,2),'-r');
plot(t,in_a(:,2),'--k');
% legend('TF','TF est','model','ref');
legend('TF est','model','ref');
hold off;

figure;
subplot(2,1,1);
plot(t,e_a(:,1));
subplot(2,1,2);
plot(t,e_a(:,2));

% [FL.p1 FL.p2]
% max(max(abs(e_a(50:end,:))))

figure; hold on;
plot(ref(1,:),ref(2,:),'--k');
plot(out(1,:),out(2,:),'-r');
hold off;
% 
% e1 = ref(1,:)-out(1,:);
% e2 = ref(2,:)-out(2,:);
% figure;
% subplot(2,1,1);
% plot(t,e1);
% subplot(2,1,2);
% plot(t,e2);

% %% Uncertain model
% 
% tolI = 0.20; % 25% accurate (inertia parameters)
% toll = 0.05; % 5% accurate (bar lengths)
% p = 4;
% 
% % Positive tolerances
% 
% uncparam_pos.g = g;
% uncparam_pos.m1 = round(m1*(1+tolI),p);
% uncparam_pos.m2 = round(m2*(1+tolI),p);
% uncparam_pos.l0 = round(l0*(1+toll),p);
% uncparam_pos.l1 = round(l1*(1+toll),p);
% uncparam_pos.l2 = round(l2*(1+toll),p);
% uncparam_pos.lg1 = round(lg1*(1+tolI),p);
% uncparam_pos.lg2 = round(lg2*(1+tolI),p);
% uncparam_pos.Jz1 = round(Jz1*(1+tolI),p);
% uncparam_pos.Jz2 = round(Jz2*(1+tolI),p);
% 
% [q_unc_pos,~,~] = ParallelRobDynamics(x0,t,param,uncparam_pos,FL, ...
%     ref,cord,T,sat,FF,sim_type);
% 
% in = ref(1:6,:);
% out = q_unc_pos;
% 
% in_a = zeros(2,length(t));
% out_a = zeros(2,length(t));
% 
% for i = 1:length(t)
%     in_a(:,i) = Qa'*in(:,i);
%     out_a(:,i) = Qa'*out(:,i);
% end
% in_a = in_a';
% out_a = out_a';
% x0_a = Qa'*x0(1:6);
% 
% data1 = iddata(out_a(:,1),in_a(:,1),T);
% data2 = iddata(out_a(:,2),in_a(:,2),T);
% 
% sys1_unc_pos = tfest(data1,2,0);
% sys2_unc_pos = tfest(data2,2,0);
% 
% G_est_unc_pos = [sys1_unc_pos 0;0 sys2_unc_pos];
% G_est_unc_pos_ss = ss(G_est_unc_pos);
% [y_est_unc_pos,~,~] = lsim(G_est_unc_pos_ss,in_a,t,[x0_a(1);x0(9); ...
%     x0_a(2);x0(11)]);
% 
% figure; hold on;
% % plot(t,y(:,1),'-b');
% plot(t,y_est_unc_pos(:,1),'-k');
% plot(t,out_a(:,1),'-r');
% plot(t,in_a(:,1),'--k');
% % legend('TF','TF est','model','ref');
% legend('TF est','model','ref');
% hold off;
% 
% figure; hold on;
% % plot(t,y(:,2),'-b');
% plot(t,y_est_unc_pos(:,2),'-k');
% plot(t,out_a(:,2),'-r');
% plot(t,in_a(:,2),'--k');
% % legend('TF','TF est','model','ref');
% legend('TF est','model','ref');
% hold off;
% 
% % Negative tolerances
% 
% uncparam_neg.g = g;
% uncparam_neg.m1 = round(m1*(1-tolI),p);
% uncparam_neg.m2 = round(m2*(1-tolI),p);
% uncparam_neg.l0 = round(l0*(1-toll),p);
% uncparam_neg.l1 = round(l1*(1-toll),p);
% uncparam_neg.l2 = round(l2*(1-toll),p);
% uncparam_neg.lg1 = round(lg1*(1-tolI),p);
% uncparam_neg.lg2 = round(lg2*(1-tolI),p);
% uncparam_neg.Jz1 = round(Jz1*(1-tolI),p);
% uncparam_neg.Jz2 = round(Jz2*(1-tolI),p);
% 
% [q_unc_neg,~,~] = ParallelRobDynamics(x0,t,param,uncparam_neg,FL, ...
%     ref,cord,T,sat,FF,sim_type);
% 
% in = ref(1:6,:);
% out = q_unc_neg;
% 
% in_a = zeros(2,length(t));
% out_a = zeros(2,length(t));
% 
% for i = 1:length(t)
%     in_a(:,i) = Qa'*in(:,i);
%     out_a(:,i) = Qa'*out(:,i);
% end
% in_a = in_a';
% out_a = out_a';
% x0_a = Qa'*x0(1:6);
% 
% data1 = iddata(out_a(:,1),in_a(:,1),T);
% data2 = iddata(out_a(:,2),in_a(:,2),T);
% 
% sys1_unc_neg = tfest(data1,2,0);
% sys2_unc_neg = tfest(data2,2,0);
% 
% G_est_unc_neg = [sys1_unc_neg 0;0 sys2_unc_neg];
% G_est_unc_neg_ss = ss(G_est_unc_neg);
% [y_est_unc_neg,~,~] = lsim(G_est_unc_neg_ss,in_a,t,[x0_a(1);x0(9); ...
%     x0_a(2);x0(11)]);
% 
% figure; hold on;
% % plot(t,y(:,1),'-b');
% plot(t,y_est_unc_neg(:,1),'-k');
% plot(t,out_a(:,1),'-r');
% plot(t,in_a(:,1),'--k');
% % legend('TF','TF est','model','ref');
% legend('TF est','model','ref');
% hold off;
% 
% figure; hold on;
% % plot(t,y(:,2),'-b');
% plot(t,y_est_unc_neg(:,2),'-k');
% plot(t,out_a(:,2),'-r');
% plot(t,in_a(:,2),'--k');
% % legend('TF','TF est','model','ref');
% legend('TF est','model','ref');
% hold off;
% 
% %% Singular values
% 
% figure;
% sigma(G_est_nom,G_est_unc_pos,G_est_unc_neg);
% legend('Nominal','Uncertainties with positive tolerances', ...
%     'Uncertainties with negative tolerances');

%------------- END OF CODE --------------