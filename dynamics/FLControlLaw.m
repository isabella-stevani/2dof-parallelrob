function [u] = FLControlLaw(x,ref,cord,param,FL,u_K,Qa,Qp,sat,FF, ...
    pos,sim_type)
% Computes robust control law.
% Inputs:
%   x: state vector
%   ref: reference signal
%   cord: actuated coordinates
%   param: model parameters
%   FL: FL control parameters
%   u_K: robust control signal
%   Qa,Qp: actuated and passive coordinates
%   sat: actuators' saturation
%   FF: feedforward enable
%   pos: time vector position
%   sim_type: simulation type
% Outputs:
%   u: control signal [2x1]

    % Reference signal coordinates
    if strcmp(cord,'xy')
        r = [ref(1);ref(2)];
        dr = [ref(7);ref(8)];
        d2r = [ref(13);ref(14)];
    elseif strcmp(cord,'theta')
        r = [ref(3);ref(5)];
        dr = [ref(9);ref(11)];
        d2r = [ref(15);ref(17)];
    end

    q = x(1:6);
    dq = x(7:12);
    
    qa = Qa'*q;
    dqa = Qa'*dq;
    ea = r-qa;
    dea = dr-dqa;
    
    if FF.on
        % Dynamics' matrices in terms of the active coordinates (previously
        % computed)
        Ma = FF.Ma(:,:,pos);
        Va = FF.Va(:,pos);
        Ga = FF.Ga(:,pos);
    else
        % Dynamic model
        [~,~,~,C,dC] = ParallelRobKinMatrix(q,dq,param,Qa,Qp);
        [M,V,G] = ParallelRobDynMatrix(q,dq,param);

        % Dynamics' matrices in terms of the active coordinates
        Ma = C'*M*C;
        Va = C'*(M*dC*dqa+V);
        Ga = C'*G;
    end
    
    if strcmp(sim_type,'default')
%         uhinf = d2r + FL.K1*dr + FL.K2*r;
%         u = sat*tanh((Va+Ga+Ma*(-FL.K1*dqa-FL.K2*qa+uhinf))/sat); %control signal
        u = sat*tanh((Va+Ga+Ma*(-FL.K1*dqa+FL.K2*u_K))/sat); %control signal
%         uhinf = r-qa;
%         u = Va+Ga+Ma*(-40*dqa+(40*2*pi*1/0.07)*r); %control signal
%         uhinf = r-qa;
%         u = Va+Ga+Ma*(-FL.K1*dqa+FL.K2*uhinf); %control signal
    elseif strcmp(sim_type,'design')
%         uhinf = d2r + FL.K1*dr+ FL.K2*r;
%         u = Va+Ga+Ma*(-FL.K1*dqa-FL.K2*qa+uhinf); %control signal
        u = Va+Ga+Ma*(-FL.K1*dqa+FL.K2*u_K); %control signal
    end

end