function [u] = RobustControlLaw(x,ref,cord,param,lambda,Qa,Qp,sat,FF,pos)
% Computes robust control law.
% Inputs:
%   x: state vector
%   ref: reference signal
%   cord: actuated coordinates
%   param: model parameters
%   lambda: FL control parameter
%   Qa,Qp: actuated and passive coordinates
%   sat: actuators' saturation
%   FF: feedforward enable
%   pos: time vector position
% Outputs:
%   u: control signal [2x1]

    % Reference signal coordinates
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
    
    qa = Qa'*q;
    dqa = Qa'*dq;
    ea = r-qa;
    dea = dr-dqa;
    
    if FF.on
        % Dynamics' matrices in terms of the active coordinates (previously
        % computed)
        Ma = FF.Ma(:,:,pos);
        Va = FF.Va(:,pos);
        Ga = FF.Ga(:,pos);
    else
        % Dynamic model
        [~,~,~,C,dC] = ParallelRobKinMatrix(q,dq,param,Qa,Qp);
        [M,V,G] = ParallelRobDynMatrix(q,dq,param);

        % Dynamics' matrices in terms of the active coordinates
        Ma = C'*M*C;
        Va = C'*(M*dC*dqa+V);
        Ga = C'*G;
    end
   
    u = sat*tanh((Va+Ga+Ma*(d2r+2*lambda*dea+lambda^2*ea))/sat); %control signal

end