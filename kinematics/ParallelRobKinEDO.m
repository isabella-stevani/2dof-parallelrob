function [dx] = ParallelRobKinEDO(t,x,param,inputfunc,lambda)
% Parallel mechanism kinematic model for Runge-Kutta method.
% Inputs:
%   t: time vector
%   x: state vector
%   param: model parameters
%   inputfunc: type of reference signal (usually constant)
%   lambda: FL control parameter
%   cord: actuated coordinates
% Outputs:
%   dx: state vector derivative [12x1]

    [r,dr,d2r] = inputfunc(t); %end-effector position reference
    
    q = x(1:6);
    dq = x(7:12);
    q0 = q(1:2);
    dq0 = dq(1:2);

    e = r-q0;
    de = dr-dq0;
    
    % Splitting model in active and passive coordinates
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
    
    % Kinematics' matrices
    [qbar,A,dA,~,~] = ParallelRobKinMatrix(q,dq,param,Qa,Qp);
    
    % Kinematic simulation model
    lambdabar = 10*lambda;
    Z = [Qa';A];
    z = [d2r + 2*lambda*de + lambda^2*e;
    -dA*dq-2*lambdabar*A*dq-lambdabar^2*qbar];
    d2q = Z\z;
    
    dx = [dq;d2q];

end