function [q] = ParallelRobKinematics(r,q12_0,param,inputfunc,lambda)
% Computes parallel mechanism kinematics.
% Inputs:
%   r: end-effector reference coordinates
%   q12_0: joints initial conditions
%   param: model parameters
%   inputfunc: type of reference signal (usually constant)
%   lambda: FL control parameter
% Outputs:
%   q: mechanism coordinates [6x1]

    
    t = [0 1]; %simulation time
    x0 = [r;q12_0;zeros(6,1)]; %initial state vector
    
    %parallel mechanism kinematic model
    [~,x] = ode45(@(t,x) ParallelRobKinEDO(t,x,param,inputfunc,lambda),t,x0);

    x = x';
    q0 = r;
    q1 = [x(3,end); x(4,end)];
    q2 = [x(5,end); x(6,end)];

    q = [q0;q1;q2]; %resultant mechanism coordinates

end
