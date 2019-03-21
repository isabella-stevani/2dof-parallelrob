function [dx] = ParallelRobDynEDO(t,x,param,uncparam,lambda,ref,cord,T)
% Parallel mechanism dynamic model for ODE method.
% Inputs:
%   t: time
%   x: state vector
%   param: model nominal parameters
%   uncparam: model real parameters
%   lambda: FL control parameter
%   ref: reference signal
%   cord: actuated coordinates
%   tvec: complete time vector
% Outputs:
%   dx: state vector derivative [12x1]

    pos = round(t/T)+1;
    ref = ref(:,pos);

    q = x(1:6);
    dq = x(7:12);
    
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
    
    lambdabar = 10*lambda;
    
    % Kinematics' matrices
    [qbarunc,Aunc,dAunc,Cunc,~] = ParallelRobKinMatrix(q,dq,uncparam,Qa,Qp);
    % Dynamics' matrices
    [Munc,Vunc,Gunc] = ParallelRobDynMatrix(q,dq,uncparam);
    % Robust control law
    u = RobustControlLaw(x,ref,param,lambda,cord);
    % Dynamic simulation model
    Z = [Cunc'*Munc;Aunc];
    z = [u+Cunc'*(-Vunc-Gunc);
    -dAunc*dq-2*lambdabar*Aunc*dq-lambdabar^2*qbarunc];
    d2q = Z\z;
    
    dx = [dq;d2q];

end