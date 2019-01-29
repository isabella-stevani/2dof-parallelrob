function [q,dq,u] = ParallelRobDynamics(x0,t,param,uncparam,lambda,ref,cord)
% Computes parallel mechanism dynamics.
% Inputs:
%   x0: initial state vector
%   t: time vector
%   param: nominal parameters for control law
%   uncparam: model real parameters
%   lambda: FL control parameter
%   ref: reference signal
%   cord: actuated coordinates
% Outputs:
%   q: mechanism coordinates [6x1]
%   dq: mechanism coordinates derivative [6x1]
%   u: control signal [2x1]
    
    T = t(2)-t(1); %sample time
    n = length(t); %number of iteractions for Runge-Kutta function
    x = ones(12,n);
    x(:,1) = x0;
    u = zeros(2,n);
    lambdabar = 10*lambda; %quick convergence for coupling equations
    func = @(x,u) ParallelRobDynRK(x,u,uncparam,lambdabar,cord); %defining
    %parallel mechanism dynamic model
    for i = 1:n-1
        u(:,i) = RobustControlLaw(x(:,i),ref(:,i),param,lambda,cord); %control
        %signal
        [x(:,i+1)] = rk4(x(:,i),u(:,i),T,func); %state vector computation
        %through 4th order Runge-Kutta method
    end
    
    q = x(1:6,:);
    dq = x(7:12,:);

end
