% Author: Isabella Stevani, Masters Student
% Master's Dissertation Project: Robust Control of a 2-DOF Parallel
% Mechanism Combining Feedback Linearization and H-Infinity Design
% Code: H-Infinity Design
% Polytechnic School of The University of Sao Paulo, Dept. of 
% Telecommunications and Control (PTC)
% E-mail address: isabella.stevani@usp.br
% Creation: Dec 2019; Last revision: 05-Apr-2020

close all; clear; clc;

%------------- BEGIN OF CODE --------------

%% Parameters

set_env; %script to set work environment

% Simulation data

tsim = 5; %simulation time [s]
t = 0:T:tsim-T; %simulation time vector
len_t = length(t);
% len_t = 1/T; t = t(len_t:end-1);
f = (1/T)*(0:(len_t/2))/len_t; %simulation frequency vector

cord = 'theta'; %reference considering end-effector coordinates
inputfunc = @RoundInput; %round reference signal

% w_in = 1:100;
% f_in = w_in/(2*pi);
% f_in = 0.5:0.5:16;
f_in = round(0.2:0.2:16,4);
% f_in = round(1:16,4);
% f_in = round(1:4:16,4);
w_in = f_in*2*pi;
% f_in = 5;
% w_in = f_in*2*pi;

x0 = x0_ref+x0_tol; %initial conditions tolerance
sim_type = 'design'; %control design simulation
FF.on = false; %without feed-forward

fig_sigma = figure; hold on;
subplot(2,1,1); hold on; subplot(2,1,2); hold on;

K = 0;

%% Simulation

% Nominal case

mag = zeros(2,length(f_in));
TF_data = zeros(2,length(f_in));

for j = 1:length(f_in)
    
    %%% Reference signal
    ref = ParallelRobDynRef(x0_ref,t,param,inputfunc,f_in(j),FL);

    [q,~,~] = ParallelRobDynamics(x0,t,param,param,FL,K, ...
        ref,cord,T,sat,FF,sim_type);

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

    in_a = zeros(2,len_t);
    out_a = zeros(2,len_t);

    for i = 1:len_t
        in_a(:,i) = Qa'*in(:,i);
        out_a(:,i) = Qa'*out(:,i);
    end
    in_a = in_a';
    out_a = out_a';
    
    df = 1/(T*len_t);
    fft_f = round((0:len_t/2)*df,4);
    fft_in = fft(in_a);
    fft_out = fft(out_a);
    fft_TF = fft_out./fft_in;
    mag2 = abs(fft_TF);
    mag1 = mag2(1:len_t/2+1,:);
    
    mag(:,j) = mag1(fft_f==f_in(j),:);
    TF_data(:,j) = fft_TF(fft_f==f_in(j),:);

end

mag_nom = mag2db(mag);
[b1,a1] = invfreqs(TF_data(1,:),w_in,0,2);
[b2,a2] = invfreqs(TF_data(2,:),w_in,0,2);
TF_nom = [tf(b1,a1) 0; 0 tf(b2,a2)];
sv_nom = mag2db(sigma(TF_nom,w_in));
sv_error_nom = sigma((TF_nom-TF)*TF^-1,w_in);
error_nom = [norm(sv_error_nom(1,:));norm(sv_error_nom(2,:))];
TF_error_nom = [TF_nom(1,1);TF_nom(2,2)];

% Uncertain case

%%% Samples between boundaries

%%%% LHS
ns = 20; %samples
% ns = 1; %samples
np = 9; %parameters
tol_LHS = lhsdesign(ns,np);

mag_unc = cell(1,ns);
TF_unc = cell(1,ns);
sv_unc = cell(1,ns);
sv_error_unc = cell(1,ns);
error_unc = cell(1,ns);

for k = 1:ns
    
    uncparam.g = g;
    uncparam.m1 = round(m1*(1+tolI*(-1+tol_LHS(k,1)*2)),p);
    uncparam.m2 = round(m2*(1+tolI*(-1+tol_LHS(k,2)*2)),p);
    uncparam.l0 = round(l0*(1+toll*(-1+tol_LHS(k,3)*2)),p);
    uncparam.l1 = round(l1*(1+toll*(-1+tol_LHS(k,4)*2)),p);
    uncparam.l2 = round(l2*(1+toll*(-1+tol_LHS(k,5)*2)),p);
    uncparam.lg1 = round(lg1*(1+tolI*(-1+tol_LHS(k,6)*2)),p);
    uncparam.lg2 = round(lg2*(1+tolI*(-1+tol_LHS(k,7)*2)),p);
    uncparam.Jz1 = round(Jz1*(1+tolI*(-1+tol_LHS(k,8)*2)),p);
    uncparam.Jz2 = round(Jz2*(1+tolI*(-1+tol_LHS(k,9)*2)),p);
    
    mag = zeros(2,length(f_in));
    TF_data = zeros(2,length(f_in));
    
    for j = 1:length(f_in)
        
        %%% Reference signal
        ref = ParallelRobDynRef(x0_ref,t,param,inputfunc,f_in(j),FL);

        [q,~,~] = ParallelRobDynamics(x0,t,param,uncparam,FL,K,ref, ...
            cord,T,sat,FF,sim_type);

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

        in_a = zeros(2,len_t);
        out_a = zeros(2,len_t);

        for i = 1:len_t
            in_a(:,i) = Qa'*in(:,i);
            out_a(:,i) = Qa'*out(:,i);
        end
        in_a = in_a';
        out_a = out_a';

        df = 1/(T*len_t);
        fft_f = round((0:len_t/2)*df,4);
        fft_in = fft(in_a);
        fft_out = fft(out_a);
        fft_TF = fft_out./fft_in;
        mag2 = abs(fft_TF);
        mag1 = mag2(1:len_t/2+1,:);

        mag(:,j) = mag1(fft_f==f_in(j),:);
        TF_data(:,j) = fft_TF(fft_f==f_in(j),:);
        
    end
    
    mag_unc{k} = mag2db(mag);
    [b1,a1] = invfreqs(TF_data(1,:),w_in,0,2);
    [b2,a2] = invfreqs(TF_data(2,:),w_in,0,2);
    TF_unc{k} = [tf(b1,a1) 0; 0 tf(b2,a2)];
    sv_unc{k} = mag2db(sigma(TF_unc{k},w_in));
    sv_error_unc{k} = sigma((TF_unc{k}-TF)*TF^-1,w_in);
    error_unc{k} = [norm(sv_error_unc{k}(1,:));norm(sv_error_unc{k}(2,:))];
    
    figure(fig_sigma);
    subplot(2,1,1);
    plot(w_in,mag_unc{k}(1,:),'-b');
    p1 = plot(w_in,mag_unc{k}(2,:),'-b');
    subplot(2,1,2);
    plot(w_in,sv_unc{k}(1,:),'-b');
    plot(w_in,sv_unc{k}(2,:),'-b');
    
end

[max_error_unc,max_error_unc_idx] = max([error_unc{:}],[],2);
aux1 = TF_unc{max_error_unc_idx(1)};
aux2 = TF_unc{max_error_unc_idx(2)};
TF_max_error_unc = [aux1(1,1);aux2(2,2)];

%%% Tolerance boundaries: Positive

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

mag = zeros(2,length(f_in));
TF_data = zeros(2,length(f_in));

for j = 1:length(f_in)
    
    %%% Reference signal
    ref = ParallelRobDynRef(x0_ref,t,param,inputfunc,f_in(j),FL);

    [q,~,~] = ParallelRobDynamics(x0,t,param,uncparam,FL,K, ...
        ref,cord,T,sat,FF,sim_type);

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

    in_a = zeros(2,len_t);
    out_a = zeros(2,len_t);

    for i = 1:len_t
        in_a(:,i) = Qa'*in(:,i);
        out_a(:,i) = Qa'*out(:,i);
    end
    in_a = in_a';
    out_a = out_a';
    
    df = 1/(T*len_t);
    fft_f = round((0:len_t/2)*df,4);
    fft_in = fft(in_a);
    fft_out = fft(out_a);
    fft_TF = fft_out./fft_in;
    mag2 = abs(fft_TF);
    mag1 = mag2(1:len_t/2+1,:);
    
    mag(:,j) = mag1(fft_f==f_in(j),:);
    TF_data(:,j) = fft_TF(fft_f==f_in(j),:);

end

mag_unc_pos = mag2db(mag);
[b1,a1] = invfreqs(TF_data(1,:),w_in,0,2);
[b2,a2] = invfreqs(TF_data(2,:),w_in,0,2);
TF_unc_pos = [tf(b1,a1) 0; 0 tf(b2,a2)];
sv_unc_pos = mag2db(sigma(TF_unc_pos,w_in));
sv_error_unc_pos = sigma((TF_unc_pos-TF)*TF^-1,w_in);
error_unc_pos = [norm(sv_error_unc_pos(1,:));norm(sv_error_unc_pos(2,:))];
TF_error_unc_pos = [TF_unc_pos(1,1);TF_unc_pos(2,2)];

%%% Tolerance boundaries: Negative

uncparam.g = g;
uncparam.m1 = round(m1*(1-tolI),p);
uncparam.m2 = round(m2*(1-tolI),p);
uncparam.l0 = round(l0*(1-toll),p);
uncparam.l1 = round(l1*(1-toll),p);
uncparam.l2 = round(l2*(1-toll),p);
uncparam.lg1 = round(lg1*(1-tolI),p);
uncparam.lg2 = round(lg2*(1-tolI),p);
uncparam.Jz1 = round(Jz1*(1-tolI),p);
uncparam.Jz2 = round(Jz2*(1-tolI),p);

mag = zeros(2,length(f_in));
TF_data = zeros(2,length(f_in));

for j = 1:length(f_in)
    
    %%% Reference signal
    ref = ParallelRobDynRef(x0_ref,t,param,inputfunc,f_in(j),FL);

    [q,~,~] = ParallelRobDynamics(x0,t,param,uncparam,FL,K, ...
        ref,cord,T,sat,FF,sim_type);

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

    in_a = zeros(2,len_t);
    out_a = zeros(2,len_t);

    for i = 1:len_t
        in_a(:,i) = Qa'*in(:,i);
        out_a(:,i) = Qa'*out(:,i);
    end
    in_a = in_a';
    out_a = out_a';
    
    df = 1/(T*len_t);
    fft_f = round((0:len_t/2)*df,4);
    fft_in = fft(in_a);
    fft_out = fft(out_a);
    fft_TF = fft_out./fft_in;
    mag2 = abs(fft_TF);
    mag1 = mag2(1:len_t/2+1,:);
    
    mag(:,j) = mag1(fft_f==f_in(j),:);
    TF_data(:,j) = fft_TF(fft_f==f_in(j),:);

end

mag_unc_neg = mag2db(mag);
[b1,a1] = invfreqs(TF_data(1,:),w_in,0,2);
[b2,a2] = invfreqs(TF_data(2,:),w_in,0,2);
TF_unc_neg = [tf(b1,a1) 0; 0 tf(b2,a2)];
sv_unc_neg = mag2db(sigma(TF_unc_neg,w_in));
sv_error_unc_neg = sigma((TF_unc_neg-TF)*TF^-1,w_in);
error_unc_neg = [norm(sv_error_unc_neg(1,:));norm(sv_error_unc_neg(2,:))];
TF_error_unc_neg = [TF_unc_neg(1,1);TF_unc_neg(2,2)];

error = [error_nom,max_error_unc,error_unc_pos,error_unc_neg];
TF_error = [TF_error_nom,TF_max_error_unc,TF_error_unc_pos, ...
    TF_error_unc_neg];

[error_i,i] = max(error,[],1);
[error_ij,j] = max(error_i);
max_error = error(i(j),j);
max_TF_error = TF_error(i(j),j);
max_sv_error = mag2db(sigma(max_TF_error,w_in));

%% Graphics

figure(fig_sigma);
subplot(2,1,1);
plot(w_in,mag_unc_pos(1,:),'-r');
p2 = plot(w_in,mag_unc_pos(2,:),'-r');
plot(w_in,mag_unc_neg(1,:),'-r');
p3 = plot(w_in,mag_unc_neg(2,:),'-r');
plot(w_in,mag_nom(1,:),'-k');
p4 = plot(w_in,mag_nom(2,:),'-k');
axis([-Inf w_in(end) -20 0]);
title('FFT');
ylabel('Magnitude (dB)');
grid on;
set(gca, 'XScale', 'log');
legend([p1,p2,p3,p4],{'Uncertain','Positive','Negative','Nominal'});
subplot(2,1,2);
plot(w_in,sv_unc_pos(1,:),'-r',w_in,sv_unc_pos(2,:),'-r');
plot(w_in,sv_unc_neg(1,:),'-r',w_in,sv_unc_neg(2,:),'-r');
plot(w_in,sv_nom(1,:),'-k',w_in,sv_nom(2,:),'-k');
axis([-Inf w_in(end) -20 0]);
title('TF');
ylabel('Magnitude');
grid on;
set(gca,'XScale','log');

%% Robustness barriers

lm = (max_TF_error-TF(1,1))*TF(1,1)^-1;
figure; sigma(lm,w_in);
axis([-Inf w_in(end) -Inf Inf]);
legend('lm','Location','Northwest');

% Parameters
delta_r = 0.07; w_r = 15; % reference tracking
delta_d = 0.1; w_d = 20; % disturbance rejection
delta_m = 0.25; w_m = 80; % measuring error rejection

% Barriers
S = 1/lm;
R = 1/(delta_r*(1-lm));
D = 1/(delta_d*(1-lm));
M = delta_m/(1+lm);

bar_s = mag2db(sigma(S,w_in));
bar_r = mag2db(sigma(R,w_in));
bar_d = mag2db(sigma(D,w_in));
bar_m = mag2db(sigma(M,w_in));

sv_nom = mag2db(sigma(TF(1,1),w_in));

for i = 1:length(w_in)
    if w_in(i) > w_r
        bar_r(i) = 0;
    end
    if w_in(i) > w_d
        bar_d(i) = 0;
    end
    if w_in(i) < w_m
        bar_m(i) = 0;
    end
end

fig_hinf = figure; hold on;
plot(w_in,bar_s,'--c');
plot(w_in,bar_r,'--r');
plot(w_in,bar_d,'--b');
plot(w_in,bar_m,'--m');
plot(w_in,sv_nom,'-k');
set(gca,'XScale','log');
axis([-Inf w_in(end) -Inf Inf]);

% Loop-shaping
s = tf('s');
G = TF(1,1);
W1 = 475/(s+15);
W2 = 15/(s+15);
W = W1*W2;
% W1 = 400/(s+25);
% W2 = (s+10)/s;
G_W = W2*G*W1;
sv_nom_W = mag2db(sigma(G_W,w_in));
% [K,CL,gamma,info] = ncfsyn(G,W1,W2);
% [K,CL,gamma,info] = loopsyn(G,G_W);
K = W;
GK = G*K;
sv_nom_GK = mag2db(sigma(GK,w_in));

figure(fig_hinf);
plot(w_in,sv_nom_W,'--k');
plot(w_in,sv_nom_GK,'*k');

legend('Stability','Tracking','Disturbance','Measurement','Plant', ...
    'Shaped Plant','Controlled Plant');

%% Comparison



save('status.mat');

% x0_a = Qa'*x0(1:6); %actuated initial conditions

% figure;plot(t_1,in_a(:,1),'--k',t_1,out_a(:,1),'-k');


% %% Parameters
% 
% set_env; %script to set work environment
% 
% % Simulation data
% 
% tsim = 3; %simulation time [s]
% t = 0:T:tsim-T; %simulation time vector
% len_t = length(t);
% f = (1/T)*(0:(len_t/2))/len_t; %simulation frequency vector
% 
% cord = 'theta'; %reference considering end-effector coordinates
% inputfunc = @RoundInput; %round reference signal
% 
% %% Nominal case
% 
% % Simulation
% 
% %%% Reference signal
% f_in = 2; %input frequency
% ref = ParallelRobDynRef(x0,t,param,inputfunc,f_in,FL);
% 
% %%% Dynamics
% x0 = x0+x0_tol; %initial conditions tolerance
% sim_type = 'design'; %control design simulation
% FF.on = false; %without feed-forward
% 
% [q,~,~] = ParallelRobDynamics(x0,t,param,param,FL, ...
%     ref,cord,T,sat,FF,sim_type);
% 
% % Transfer function estimation
% 
% %%% Raw data
% in = ref(1:6,:);
% out = q;
% 
% %%% Actuated data
% if strcmp(cord,'xy')
%     Qa = [1 0;
%         0 1;
%         0 0;
%         0 0;
%         0 0;
%         0 0];
% elseif strcmp(cord,'theta')
%     Qa = [0 0;
%         0 0;
%         1 0;
%         0 0;
%         0 1;
%         0 0];
% end
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
% 
% x0_a = Qa'*x0(1:6); %actuated initial conditions
% 
% e_a = in_a-out_a;
% 
% % data = iddata(out_a,in_a,T);
% data1 = iddata(out_a(:,1),in_a(:,1),T); %first actuated variable data
% data2 = iddata(out_a(:,2),in_a(:,2),T); %second actuated variable data
% 
% %%% TF Estimation
% % tfest_opt = tfestOptions;
% % tfest_opt.InitializeMethod = 'all';
% % sys = tfest(data,2,0);
% sys1_nom = tfest(data1,2,0); %estimated TF for the first actuated variable
% sys2_nom = tfest(data2,2,0); %estimated TF for the second actuated variable
% 
% %%%% Estimated TF response
% G_est_nom = [sys1_nom 0;0 sys2_nom];
% G_est_nom_ss = ss(G_est_nom);
% [y_est_nom,~,~] = lsim(G_est_nom_ss,in_a,t,[x0_a(1);x0(9);x0_a(2); ...
%     x0(11)]);
% 
%%% Frequency analysis
% FFT_in_1 = fft(in_a(:,1));
% P2_in_1 = abs(FFT_in_1/len_t);
% P1_in_1 = P2_in_1(1:len_t/2+1);
% P1_in_1(2:end-1) = 2*P1_in_1(2:end-1);
% FFT_out_1 = fft(out_a(:,1));
% P2_out_1 = abs(FFT_out_1/len_t);
% P1_out_1 = P2_out_1(1:len_t/2+1);
% P1_out_1(2:end-1) = 2*P1_out_1(2:end-1);
% 
% FFT_in_2 = fft(in_a(:,2));
% P2_in_2 = abs(FFT_in_2/len_t);
% P1_in_2 = P2_in_2(1:len_t/2+1);
% P1_in_2(2:end-1) = 2*P1_in_2(2:end-1);
% FFT_out_2 = fft(out_a(:,2));
% P2_out_2 = abs(FFT_out_2/len_t);
% P1_out_2 = P2_out_2(1:len_t/2+1);
% P1_out_2(2:end-1) = 2*P1_out_2(2:end-1);
% 
% figure;
% subplot(2,1,1);
%plot(f(len_1:end),P1_in_1,'--k',f(len_1:end),P1_out_1,'-k');
% plot(fft_f,P1_in_1,'--k',fft_f,P1_out_1,'-k');
% axis([0 50 -Inf Inf]);
% subplot(2,1,2);
% plot(f,P1_in_2,'--k',f,P1_out_2,'-k');
% legend('Reference','Output');
% axis([0 50 -Inf Inf]);
% 
% % Nominal model
% 
% %%% Continuous linearized model
% s = tf('s');
% G = FL.K2/(s^2+FL.K1*s+FL.K2)*eye(2);
% G_ss = ss(G);
% [y,~,~] = lsim(G_ss,in_a,t,[x0_a(1);x0(9);x0_a(2);x0(11)]);
% 
% % Plots
% 
% % figure;
% % sigma(G,sys(1,1),sys(2,2));
% % legend on;
% % figure;
% % sigma(G,sys1);
% % legend on;
% 
% % figure; hold on;
% % % plot(t,y(:,1),'-b');
% % plot(t,y_est_nom(:,1),'-k');
% % plot(t,out_a(:,1),'-r');
% % plot(t,in_a(:,1),'--k');
% % % legend('TF','TF est','model','ref');
% % legend('TF est','model','ref');
% % hold off;
% % 
% % figure; hold on;
% % % plot(t,y(:,2),'-b');
% % plot(t,y_est_nom(:,2),'-k');
% % plot(t,out_a(:,2),'-r');
% % plot(t,in_a(:,2),'--k');
% % % legend('TF','TF est','model','ref');
% % legend('TF est','model','ref');
% % hold off;
% % 
% % figure;
% % subplot(2,1,1);
% % plot(t,e_a(:,1));
% % subplot(2,1,2);
% % plot(t,e_a(:,2));
% 
% % [FL.p1 FL.p2]
% % max(max(abs(e_a(50:end,:))))
% 
% % figure; hold on;
% % plot(ref(1,:),ref(2,:),'--k');
% % plot(out(1,:),out(2,:),'-r');
% % hold off;
% % 
% % e1 = ref(1,:)-out(1,:);
% % e2 = ref(2,:)-out(2,:);
% % figure;
% % subplot(2,1,1);
% % plot(t,e1);
% % subplot(2,1,2);
% % plot(t,e2);
% 
% %% Uncertain model
% 
% FF.on = true;
% 
% % Tolerance boundaries
% 
% %%% Positive
% 
% uncparam.g = g;
% uncparam.m1 = round(m1*(1+tolI),p);
% uncparam.m2 = round(m2*(1+tolI),p);
% uncparam.l0 = round(l0*(1+toll),p);
% uncparam.l1 = round(l1*(1+toll),p);
% uncparam.l2 = round(l2*(1+toll),p);
% uncparam.lg1 = round(lg1*(1+tolI),p);
% uncparam.lg2 = round(lg2*(1+tolI),p);
% uncparam.Jz1 = round(Jz1*(1+tolI),p);
% uncparam.Jz2 = round(Jz2*(1+tolI),p);
% 
% [q_unc_pos,~,~] = ParallelRobDynamics(x0,t,param,uncparam,FL, ...
%     ref,cord,T,sat,FF,sim_type);
% 
% %%%% TF estimation
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
% %%% Frequency analysis
% FFT_out_1_unc_pos = fft(out_a(:,1));
% P2_out_1_unc_pos = abs(FFT_out_1_unc_pos/len_t);
% P1_out_1_unc_pos = P2_out_1_unc_pos(1:len_t/2+1);
% P1_out_1_unc_pos(2:end-1) = 2*P1_out_1_unc_pos(2:end-1);
% 
% FFT_out_2_unc_pos = fft(out_a(:,2));
% P2_out_2_unc_pos = abs(FFT_out_2_unc_pos/len_t);
% P1_out_2_unc_pos = P2_out_2_unc_pos(1:len_t/2+1);
% P1_out_2_unc_pos(2:end-1) = 2*P1_out_2_unc_pos(2:end-1);
% 
% % figure; hold on;
% % plot(t,y(:,1),'-b');
% % plot(t,y_est_unc_pos(:,1),'-k');
% % plot(t,out_a(:,1),'-r');
% % plot(t,in_a(:,1),'--k');
% % legend('TF','TF est','model','ref');
% % % legend('TF est','model','ref');
% % hold off;
% % 
% % figure; hold on;
% % plot(t,y(:,2),'-b');
% % plot(t,y_est_unc_pos(:,2),'-k');
% % plot(t,out_a(:,2),'-r');
% % plot(t,in_a(:,2),'--k');
% % legend('TF','TF est','model','ref');
% % % legend('TF est','model','ref');
% % hold off;
% 
% %%% Negative
% 
% uncparam.m1 = round(m1*(1-tolI),p);
% uncparam.m2 = round(m2*(1-tolI),p);
% uncparam.l0 = round(l0*(1-toll),p);
% uncparam.l1 = round(l1*(1-toll),p);
% uncparam.l2 = round(l2*(1-toll),p);
% uncparam.lg1 = round(lg1*(1-tolI),p);
% uncparam.lg2 = round(lg2*(1-tolI),p);
% uncparam.Jz1 = round(Jz1*(1-tolI),p);
% uncparam.Jz2 = round(Jz2*(1-tolI),p);
% 
% [q_unc_neg,~,~] = ParallelRobDynamics(x0,t,param,uncparam,FL, ...
%     ref,cord,T,sat,FF,sim_type);
% 
% %%%% TF estimation
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
% %%% Frequency analysis
% FFT_out_1_unc_neg = fft(out_a(:,1));
% P2_out_1_unc_neg = abs(FFT_out_1_unc_neg/len_t);
% P1_out_1_unc_neg = P2_out_1_unc_neg(1:len_t/2+1);
% P1_out_1_unc_neg(2:end-1) = 2*P1_out_1_unc_neg(2:end-1);
% 
% FFT_out_2_unc_neg = fft(out_a(:,2));
% P2_out_2_unc_neg = abs(FFT_out_2_unc_neg/len_t);
% P1_out_2_unc_neg = P2_out_2_unc_neg(1:len_t/2+1);
% P1_out_2_unc_neg(2:end-1) = 2*P1_out_2_unc_neg(2:end-1);
% 
% % figure; hold on;
% % plot(t,y(:,1),'-b');
% % plot(t,y_est_unc_neg(:,1),'-k');
% % plot(t,out_a(:,1),'-r');
% % plot(t,in_a(:,1),'--k');
% % legend('TF','TF est','model','ref');
% % % legend('TF est','model','ref');
% % hold off;
% % 
% % figure; hold on;
% % plot(t,y(:,2),'-b');
% % plot(t,y_est_unc_neg(:,2),'-k');
% % plot(t,out_a(:,2),'-r');
% % plot(t,in_a(:,2),'--k');
% % legend('TF','TF est','model','ref');
% % % legend('TF est','model','ref');
% % hold off;
% 
% % Samples between boundaries
% 
% %%% LHS
% ns = 2; %samples
% np = 9; %parameters
% tol_LHS = lhsdesign(ns,np);
% 
% fig_sv = figure; hold on;
% fig_fft = figure; hold on;
% subplot(2,1,1); hold on; subplot(2,1,2); hold on;
% for i = 1:ns
%     uncparam.m1 = round(m1*(1+tolI*(-1+tol_LHS(i,1)*2)),p);
%     uncparam.m2 = round(m2*(1+tolI*(-1+tol_LHS(i,2)*2)),p);
%     uncparam.l0 = round(l0*(1+toll*(-1+tol_LHS(i,3)*2)),p);
%     uncparam.l1 = round(l1*(1+toll*(-1+tol_LHS(i,4)*2)),p);
%     uncparam.l2 = round(l2*(1+toll*(-1+tol_LHS(i,5)*2)),p);
%     uncparam.lg1 = round(lg1*(1+tolI*(-1+tol_LHS(i,6)*2)),p);
%     uncparam.lg2 = round(lg2*(1+tolI*(-1+tol_LHS(i,7)*2)),p);
%     uncparam.Jz1 = round(Jz1*(1+tolI*(-1+tol_LHS(i,8)*2)),p);
%     uncparam.Jz2 = round(Jz2*(1+tolI*(-1+tol_LHS(i,9)*2)),p);
%     [q_unc,~,~] = ParallelRobDynamics(x0,t,param,uncparam,FL,ref,cord, ...
%         T,sat,FF,sim_type);
%                                   
%     in = ref(1:6,:);
%     out = q_unc;
% 
%     in_a = zeros(2,length(t));
%     out_a = zeros(2,length(t));
% 
%     for i = 1:length(t)
%         in_a(:,i) = Qa'*in(:,i);
%         out_a(:,i) = Qa'*out(:,i);
%     end
%     in_a = in_a';
%     out_a = out_a';
%     x0_a = Qa'*x0(1:6);
% 
%     data1 = iddata(out_a(:,1),in_a(:,1),T);
%     data2 = iddata(out_a(:,2),in_a(:,2),T);
% 
%     sys1_unc = tfest(data1,2,0);
%     sys2_unc = tfest(data2,2,0);
% 
%     G_est_unc = [sys1_unc 0;0 sys2_unc];
%     G_est_unc_ss = ss(G_est_unc);
% 
%     figure(fig_sv); sigma(G_est_unc,'-b');
%     
%     %%% Frequency analysis
%     FFT_out_1 = fft(out_a(:,1));
%     P2_out_1 = abs(FFT_out_1/len_t);
%     P1_out_1 = P2_out_1(1:len_t/2+1);
%     P1_out_1(2:end-1) = 2*P1_out_1(2:end-1);
% 
%     FFT_out_2 = fft(out_a(:,2));
%     P2_out_2 = abs(FFT_out_2/len_t);
%     P1_out_2 = P2_out_2(1:len_t/2+1);
%     P1_out_2(2:end-1) = 2*P1_out_2(2:end-1);
% 
%     figure(fig_fft);
%     subplot(2,1,1);
%     plot(f,P1_out_1,'-b');
%     subplot(2,1,2);
%     plot(f,P1_out_2,'-b');
% end
% 
% figure(fig_sv);
% sigma(G,'-k',G_est_unc_pos,'-r',G_est_unc_neg,'-r');
% hold off;
% 
% figure(fig_fft);
% subplot(2,1,1);
% plot(f,P1_out_1_unc_pos,'-r');
% plot(f,P1_out_1_unc_neg,'-r');
% plot(f,P1_in_1,'--k');
% axis([0 50 -Inf Inf]);
% subplot(2,1,2);
% plot(f,P1_out_2_unc_pos,'-r');
% plot(f,P1_out_2_unc_neg,'-r');
% plot(f,P1_in_2,'--k');
% axis([0 50 -Inf Inf]);
% 
% % %% Singular values
% % 
% % figure;
% % sigma(G,G_est_unc_pos,G_est_unc_neg);
% % legend('Nominal','Uncertainties with positive tolerances', ...
% %     'Uncertainties with negative tolerances');

%------------- END OF CODE --------------