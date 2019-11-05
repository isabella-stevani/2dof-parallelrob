function [ref] = ParallelRobDynRef(x0,t,param,inputfunc,lambda)
% Reference signal computation.
% Inputs:
%   x0: initial state vector
%   t: time vector (not used)
%   param: nominal parameters for control law
%   inputfunc: type of reference signal
%   lambda: FL control parameter
%   cord: actuated coordinates
% Outputs:
%   ref: reference signal [2x1]
    
    % Reference state vector through kinematics
    [~,x] = ode45(@(t,x) ParallelRobKinEDO(t,x,param,inputfunc,lambda),t,x0);
    x = x'; t = t';
    % State vector derivative
    dx = zeros(size(x));
    for i = 1:length(t)
        dx(:,i) = ParallelRobKinEDO(t(i),x(:,i),param,inputfunc,lambda);
    end

    r = x(1:6,:);
    dr = x(7:12,:);
    d2r = dx(7:12,:);
    
    ref = [r;dr;d2r];
    
end