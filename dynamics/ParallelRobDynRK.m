function [dx] = ParallelRobDynRK(x,u,param,lambdabar,cord)
% Parallel mechanism dynamic model for Runge-Kutta method.
% Inputs:
%   x: state vector
%   u: control signal
%   param: model parameters
%   lambdabar: parameter for convergence of the coupling equations
%   cord: actuated coordinates
% Outputs:
%   dx: state vector derivative [12x1]

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
    
    % Kinematics' matrices
    [qbarunc,Aunc,dAunc,Cunc,~] = ParallelRobKinMatrix(q,dq,param,Qa,Qp);
    % Dynamics' matrices
    [Munc,Vunc,Gunc] = ParallelRobDynMatrix(q,dq,param);
    % Dynamic simulation model
    Z = [Cunc'*Munc;Aunc];
    z = [u+Cunc'*(-Vunc-Gunc);
    -dAunc*dq-2*lambdabar*Aunc*dq-lambdabar^2*qbarunc];
    d2q = Z\z;
   
    dx = [dq;d2q];

end