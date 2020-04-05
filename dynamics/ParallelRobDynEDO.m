function [out] = ParallelRobDynEDO(t,x,param,uncparam,FL,K,ref,Qa,Qp, ...
    cord,T,sat,outvar,FF,sim_type)
% Parallel mechanism dynamic model for ODE method.
% Inputs:
%   t: time
%   x: state vector
%   param: model nominal parameters
%   uncparam: model real parameters
%   FL: FL control parameters
%   K: robust control parameters
%   ref: reference signal
%   Qa,Qp: actuated and passive coordinates
%   cord: actuated coordinates
%   T: sample time
%   sat: actuators' saturation
%   outvar: output variable, 'x' or 'u'
%   FF: feedforward enable
%   sim_type: simulation type
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
    
    %%% Reference signal coordinates
    if strcmp(cord,'xy')
        r = [ref(1);ref(2)];
        dr = [ref(7);ref(8)];
        d2r = [ref(13);ref(14)];
    elseif strcmp(cord,'theta')
        r = [ref(3);ref(5)];
        dr = [ref(9);ref(11)];
        d2r = [ref(15);ref(17)];
    end

    q = x(1:6);
    dq = x(7:12);
    u_K = x(13:14);
    du_K = x(15:16);
    
    qa = Qa'*q;
    dqa = Qa'*dq;
    ea = r-qa;
    dea = dr-dqa;
    
    %%% Robust control signal
    num_K = K.num{:};
    den_K = K.den{:};
    if length(den_K) == 1
        u_K = ea;
        du_K = [0;0];
        d2u_K = [0;0];
    else
        n2 = num_K(1);
        n1 = num_K(2);
        n0 = num_K(3);
        d2 = den_K(1);
        d1 = den_K(2);
        d0 = den_K(3);
        d2u_K = (1/d2)*(n1*dea+n0*ea-d1*du_K-d0*u_K);
%         d2u_K = n1*dea+n0*ea-d1*du_K-d0*u_K;
    end
    
    % FL control law
    u = FLControlLaw(x,ref,cord,param,FL,u_K,Qa,Qp,sat,FF,pos, ...
        sim_type);
    % Dynamic simulation model
    Z = [Cunc'*Munc;Aunc];
    z = [u+Cunc'*(-Vunc-Gunc);
    -dAunc*dq-2*FL.lambdabar*Aunc*dq-FL.lambdabar^2*qbarunc];
    d2q = Z\z;
    dx = [dq;d2q;du_K;d2u_K];
    
    if strcmp(outvar,'x')
        out = dx;
    elseif strcmp(outvar,'u')
        out = u;
    end
end