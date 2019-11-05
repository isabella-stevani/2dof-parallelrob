function [u] = RobustControlLaw(x,ref,cord,param,lambda,Qa,Qp,sat,FF)
% Computes robust control law.
% Inputs:
%   x: state vector
%   ref: reference signal
%   cord: actuated coordinates
%   param: model parameters
%   lambda: FL control parameter
%   Qa,Qp: actuated and passive coordinates
%   sat: actuators' saturation
%   FF: feedforward enable
% Outputs:
%   u: control signal [2x1]

    % Reference signal coordinates
    if strcmp(cord,'xy')
        r = [ref(1,:);ref(2,:)];
        dr = [ref(7,:);ref(8,:)];
        d2r = [ref(13,:);ref(14,:)];
    elseif strcmp(cord,'theta')
        r = [ref(3,:);ref(5,:)];
        dr = [ref(9,:);ref(11,:)];
        d2r = [ref(15,:);ref(17,:)];
    end

    q = x(1:6);
    dq = x(7:12);
    
%     r = ref(1:2);
%     dr = ref(3:4);
%     d2r = ref(5:6);
    
    qa = Qa'*q;
    dqa = Qa'*dq;
    ea = r-qa;
    dea = dr-dqa;
    
    if FF.on
        q = ref(1:6);
        dq = ref(7:12);
        dqa = Qa'*dq;
    end
    
    % Dynamic model
    [~,~,~,C,dC] = ParallelRobKinMatrix(q,dq,param,Qa,Qp);
    [M,V,G] = ParallelRobDynMatrix(q,dq,param);
    
    % Dynamics' matrices in terms of the active coordinates
    Ma = C'*M*C;
    Va = C'*(M*dC*dqa+V);
    Ga = C'*G;
    u = sat*tanh((Va+Ga+Ma*(d2r+2*lambda*dea+lambda^2*ea))/sat); %control signal

end