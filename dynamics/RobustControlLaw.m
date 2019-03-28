function [u] = RobustControlLaw(x,ref,param,lambda,cord,sat)
% Computes robust control law.
% Inputs:
%   x: state vector
%   ref: reference signal
%   param: model parameters
%   lambda: FL control parameter
%   cord: actuated coordinates
%   sat: actuators' saturation
% Outputs:
%   u: control signal [2x1]

    q = x(1:6);
    dq = x(7:12);
    
    r = ref(1:2);
    dr = ref(3:4);
    d2r = ref(5:6);
    
    % Splitting model in active and passive coordinates
    if strcmp(cord,'xy')
        Qa = [1 0;
              0 1;
              0 0;
              0 0;
              0 0;
              0 0];

        Qp = [0 0 0 0;
               0 0 0 0;
               1 0 0 0;
               0 1 0 0;
               0 0 1 0;
               0 0 0 1];
    elseif strcmp(cord,'theta')
        Qa = [0 0;
              0 0;
              1 0;
              0 0;
              0 1;
              0 0];

        Qp = [1 0 0 0;
               0 1 0 0;
               0 0 0 0;
               0 0 1 0;
               0 0 0 0;
               0 0 0 1]; 
    end
    
    qa = Qa'*q;
    dqa = Qa'*dq;
    ea = r-qa;
    dea = dr-dqa;
    
    % Dynamic model
    [~,~,~,C,dC] = ParallelRobKinMatrix(q,dq,param,Qa,Qp);
    [M,V,G] = ParallelRobDynMatrix(q,dq,param);
    
    % Dynamics' matrices in terms of the active coordinates
    Ma = C'*M*C;
    Va = C'*(M*dC*dqa+V);
    Ga = C'*G;
    u = sat*tanh((Va+Ga+Ma*(d2r+2*lambda*dea+lambda^2*ea))/sat); %control signal

end