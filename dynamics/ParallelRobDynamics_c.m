function [q,dq,u] = ParallelRobDynamics_c(x0,t,param,uncparam,lambda, ...
    ref,cord,T,sat)
% Computes parallel mechanism continuous dynamics.
% Inputs:
%   x0: initial state vector
%   t: time vector
%   param: nominal parameters for control law
%   uncparam: model real parameters
%   lambda: FL control parameter
%   ref: reference signal
%   cord: actuated coordinates
%   T: sample time
%   sat: actuators' saturation
% Outputs:
%   q: mechanism coordinates [6x1]
%   dq: mechanism coordinates derivative [6x1]
%   u: control signal [2x1]
    
    lambdabar = 10*lambda; %quick convergence for coupling equations
    [~,x] = ode45(@(t,x) ParallelRobDynEDO(t,x,param,uncparam,lambda, ...
        ref,cord,T,sat,'x'),t,x0);
    
    x = x';
    q = x(1:6,:);
    dq = x(7:12,:);
    u = zeros(2,length(t));
    
    for i = 1:length(t)
        u(:,i) = ParallelRobDynEDO(t(i),x(:,i),param,uncparam,lambda, ...
        ref,cord,T,sat,'u');
    end

end
