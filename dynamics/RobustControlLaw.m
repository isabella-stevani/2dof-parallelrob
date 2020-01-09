function [u] = RobustControlLaw(x,ref,cord,param,FL,Qa,Qp,sat,FF, ...
    pos,sim_type)
% Computes robust control law.
% Inputs:
%   x: state vector
%   ref: reference signal
%   cord: actuated coordinates
%   param: model parameters
%   FL: FL control parameters
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
    
    K1 = -(FL.p1+FL.p2);
    K2 = FL.p1*FL.p2;
    if strcmp(sim_type,'default')
        uhinf = d2r + 2*FL*dr+ FL^2*r;
        u = sat*tanh((Va+Ga+Ma*(-K1*dqa-K2*qa+uhinf))/sat); %control signal
    elseif strcmp(sim_type,'design')
%         uhinf = d2r + 2*lambda*dr+ lambda^2*r;
%         u = Va+Ga+Ma*(-2*lambda*dqa-lambda^2*qa+uhinf); %control signal
%         u = Va+Ga+Ma*(-2*lambda*dqa-lambda^2*qa)+r; %control signal
%         u = Va+Ga+Ma*(-2*lambda*dqa-lambda^2*qa+r); %control signal
%         uhinf = d2r + 2*lambda*dr+ lambda^2*r;
%         u = Va+Ga+Ma*(-2*lambda*dqa-lambda^2*qa+uhinf); %control signal
        uhinf = r-qa;
        u = Va+Ga+Ma*(-K1*dqa+K2*uhinf); %control signal
    end

end