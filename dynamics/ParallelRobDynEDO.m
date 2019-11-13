function [out] = ParallelRobDynEDO(t,x,param,uncparam,lambda,lambdabar, ...
    ref,Qa,Qp,cord,T,sat,outvar,FF)
% Parallel mechanism dynamic model for ODE method.
% Inputs:
%   t: time
%   x: state vector
%   param: model nominal parameters
%   uncparam: model real parameters
%   lambda: FL control parameter
%   lambdabar: FL model convergence parameter
%   ref: reference signal
%   Qa,Qp: actuated and passive coordinates
%   cord: actuated coordinates
%   T: sample time
%   sat: actuators' saturation
%   outvar: output variable, 'x' or 'u'
% Outputs:
%   out could be
%   dx: state vector derivative [12x1]
%   or
%   u: control signal [2x1]

    pos = round(t/T)+1;
    ref = ref(:,pos);

    q = x(1:6);
    dq = x(7:12);
    
    % Kinematics' matrices
    [qbarunc,Aunc,dAunc,Cunc,~] = ParallelRobKinMatrix(q,dq,uncparam,Qa,Qp);
    % Dynamics' matrices
    [Munc,Vunc,Gunc] = ParallelRobDynMatrix(q,dq,uncparam);
    % Robust control law
    u = RobustControlLaw(x,ref,cord,param,lambda,Qa,Qp,sat,FF,pos);
    % Dynamic simulation model
    Z = [Cunc'*Munc;Aunc];
    z = [u+Cunc'*(-Vunc-Gunc);
    -dAunc*dq-2*lambdabar*Aunc*dq-lambdabar^2*qbarunc];
    d2q = Z\z;
    dx = [dq;d2q];
    
    if strcmp(outvar,'x')
        out = dx;
    elseif strcmp(outvar,'u')
        out = u;
    end
end