% Author: Isabella Stevani, Masters Student
% Master's Dissertation Project: Robust Control of a 2-DOF Parallel
% Mechanism Combining Feedback Linearization and H-Infinity Design
% Code: Continuous Control Design and Simulation
% Polytechnic School of The University of Sao Paulo, Dept. of 
% Telecommunications and Control (PTC)
% E-mail address: isabella.stevani@usp.br
% Creation: Aug 2018; Last revision: 11-Feb-2020

close all; clear; clc;

%------------- BEGIN OF CODE --------------

%% Parameters

set_env; %script to set work environment

% Simulation data

tsim = 3; %simulation time [s]
t = 0:T:tsim-T; %simulation time vector
len_t = length(t);
f = (1/T)*(0:(len_t/2))/len_t; %simulation frequency vector
f_in = 2; %input frequency

cord = 'theta'; %reference considering end-effector coordinates
inputfunc = @RoundInput; %round reference signal

%%% Reference signal
ref = ParallelRobDynRef(x0_ref,t,param,inputfunc,f_in,FL);

x0 = x0_ref+x0_tol; %initial conditions tolerance
sim_type = 'default'; %default simulation

K = 0;

%% Nominal simulation

% Without feedforward
%%% Dynamics
FF.on = false;
[q,dq,u] = ParallelRobDynamics(x0,t,param,param,FL,K,ref,cord,T, ...
    sat,FF,sim_type);

% With feedforward
%%% Dynamics
FF.on = true;
[q_FF,dq_FF,u_FF] = ParallelRobDynamics(x0,t,param,param,FL,K, ...
    ref,cord,T,sat,FF,sim_type);

% States
figure; leg = {'Without FF','With FF'};
ylab = {'x [m]','y [m]','t11 [rad]','t12 [rad]','t21 [rad]','t22 [rad]'};
xlab = 'Time [s]';
tlab = 'States';
for i = 1:6
    subplot(3,2,i); hold on;
    plot(t,q(i,:),'-b');
    plot(t,q_FF(i,:),'-r');
    hold off;
    grid on;
    xl = xlabel(xlab);
    yl = ylabel(ylab{i});
    set(xl,'FontSize',16);
    set(yl,'FontSize',16);
end
lgd = legend(leg);
set(lgd,'FontSize',12);
supertitle(tlab,18);
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 2;
end
set(gcf,'PaperPositionMode','auto');
saveas(gcf,['images/continuous/cntr_nomsim_states','.eps'], 'epsc');
saveas(gcf,'images/continuous/cntr_nomsim_states.png');
% 
% Control signal
figure; leg = {'Without FF','With FF'};
xlab = {'','Time [s]'};
ylab = {'u1 [N.m]','u2 [N.m]'};
tlab = 'Control Signals';
for i = 1:2
    subplot(2,1,i); hold on;
    plot(t,u(i,:),'-b');
    plot(t,u_FF(i,:),'-r');
    grid on;
    hold off;
    xl = xlabel(xlab{i});
    yl = ylabel(ylab{i});
    set(xl,'FontSize',16);
    set(yl,'FontSize',16);
end
lgd = legend(leg);
set(lgd,'FontSize',12);
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 2;
end
supertitle(tlab,18);
set(gcf,'PaperPositionMode','auto');
saveas(gcf,['images/continuous/cntr_nomsim_cntrsig','.eps'], 'epsc');
saveas(gcf,'images/continuous/cntr_nomsim_cntrsig.png');

% Actuated states errors
if strcmp(cord,'xy')
    r = [ref(1,:);ref(2,:)];
    e = r-[q(1,:);q(2,:)];
    e_FF = r-[q_FF(1,:);q_FF(2,:)];
elseif strcmp(cord,'theta')
    r = [ref(3,:);ref(5,:)];
    e = r-[q(3,:);q(5,:)];
    e_FF = r-[q_FF(3,:);q_FF(5,:)];
end
figure; leg = {'Without FF','With FF'};
xlab = {'','Time [s]'};
ylab = {'e1 [rad]','e2 [rad]'};
tlab = 'Actuated Angle Errors';
for i = 1:2
    subplot(2,1,i); hold on;
    plot(t,e(i,:),'-b');
    plot(t,e_FF(i,:),'-r');
    grid on;
    hold off;
    xl = xlabel(xlab{i});
    yl = ylabel(ylab{i});
    set(xl,'FontSize',16);
    set(yl,'FontSize',16);
end
lgd = legend(leg);
set(lgd,'FontSize',12);
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 2;
end
supertitle(tlab,18);
set(gcf,'PaperPositionMode','auto');
saveas(gcf,['images/continuous/cntr_nomsim_error','.eps'], 'epsc');
saveas(gcf,'images/continuous/cntr_nomsim_error.png');

%% Uncertain simulation

n = 10; %samples

fig_x = figure;
fig_u = figure;
fig_e = figure;

for i = 1:n
    
    uncparam.g = g;
    uncparam.m1 = round(m1*(1+tolI*(-1+rand(1)*2)),p);
    uncparam.m2 = round(m2*(1+tolI*(-1+rand(1)*2)),p);
    uncparam.l0 = round(l0*(1+toll*(-1+rand(1)*2)),p);
    uncparam.l1 = round(l1*(1+toll*(-1+rand(1)*2)),p);
    uncparam.l2 = round(l2*(1+toll*(-1+rand(1)*2)),p);
    uncparam.lg1 = round(lg1*(1+tolI*(-1+rand(1)*2)),p);
    uncparam.lg2 = round(lg2*(1+tolI*(-1+rand(1)*2)),p);
    uncparam.Jz1 = round(Jz1*(1+tolI*(-1+rand(1)*2)),p);
    uncparam.Jz2 = round(Jz2*(1+tolI*(-1+rand(1)*2)),p);
    
    % Without feedforward
    %%% Dynamics
    FF.on = false;
    [q,dq,u] = ParallelRobDynamics(x0,t,param,uncparam,FL,K,ref,cord, ...
        T,sat,FF,sim_type);

    % With feedforward
    %%% Dynamics
    FF.on = true;
    [q_FF,dq_FF,u_FF] = ParallelRobDynamics(x0,t,param,uncparam,FL,K, ...
        ref,cord,T,sat,FF,sim_type);
    
    % States
    figure(fig_x);
    for i = 1:6
        subplot(3,2,i); hold on;
        plot(t,q(i,:),'-b');
        plot(t,q_FF(i,:),'-r');
    end
    
    % Control signal
    figure(fig_u);
    for i = 1:2
        subplot(2,1,i); hold on;
        plot(t,u(i,:),'-b');
        plot(t,u_FF(i,:),'-r');
    end

    % Actuated states errors
    if strcmp(cord,'xy')
        r = [ref(1,:);ref(2,:)];
        e = r-[q(1,:);q(2,:)];
        e_FF = r-[q_FF(1,:);q_FF(2,:)];
    elseif strcmp(cord,'theta')
        r = [ref(3,:);ref(5,:)];
        e = r-[q(3,:);q(5,:)];
        e_FF = r-[q_FF(3,:);q_FF(5,:)];
    end
    figure(fig_e);
    for i = 1:2
        subplot(2,1,i); hold on;
        plot(t,e(i,:),'-b');
        plot(t,e_FF(i,:),'-r');
    end
    
end

% Graphics

%%% States
leg = {'Without FF','With FF'};
ylab = {'x [m]','y [m]','t11 [rad]','t12 [rad]','t21 [rad]','t22 [rad]'};
xlab = 'Time [s]';
tlab = 'States';
figure(fig_x); hold off;
for i = 1:6
    subplot(3,2,i);
    grid on;
    xl = xlabel(xlab);
    yl = ylabel(ylab{i});
    set(xl,'FontSize',16);
    set(yl,'FontSize',16);
end
lgd = legend(leg);
set(lgd,'FontSize',12);
supertitle(tlab,18);
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
    lines(i).LineWidth = 2;
end
set(gcf,'PaperPositionMode','auto');
saveas(gcf,['images/continuous/cntr_uncsim_states','.eps'], 'epsc');
saveas(gcf,'images/continuous/cntr_uncsim_states.png');

%%% Control signal
leg = {'Without FF','With FF'};
xlab = {'','Time [s]'};
ylab = {'u1 [N.m]','u2 [N.m]'};
tlab = 'Control Signals';
figure(fig_u); hold off;
for i = 1:2
    subplot(2,1,i)
    grid on;
    xl = xlabel(xlab{i});
    yl = ylabel(ylab{i});
    set(xl,'FontSize',16);
    set(yl,'FontSize',16);
end
lgd = legend(leg);
set(lgd,'FontSize',12);
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
    lines(i).LineWidth = 2;
end
supertitle(tlab,18);
set(gcf,'PaperPositionMode','auto');
saveas(gcf,['images/continuous/cntr_uncsim_cntrsig','.eps'], 'epsc');
saveas(gcf,'images/continuous/cntr_uncsim_cntrsig.png');

%%% Actuated states errors
leg = {'Without FF','With FF'};
xlab = {'','Time [s]'};
ylab = {'e1 [rad]','e2 [rad]'};
tlab = 'Actuated Angle Errors';
figure(fig_e); hold off;
for i = 1:2
    subplot(2,1,i)
    grid on;
    xl = xlabel(xlab{i});
    yl = ylabel(ylab{i});
    set(xl,'FontSize',16);
    set(yl,'FontSize',16);
end
lgd = legend(leg);
set(lgd,'FontSize',12);
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
    lines(i).LineWidth = 2;
end
supertitle(tlab,18);
set(gcf,'PaperPositionMode','auto');
saveas(gcf,['images/continuous/cntr_uncsim_error','.eps'], 'epsc');
saveas(gcf,'images/continuous/cntr_uncsim_error.png');

%------------- END OF CODE --------------