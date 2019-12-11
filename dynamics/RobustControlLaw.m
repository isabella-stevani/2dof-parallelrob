function [u] = RobustControlLaw(x,ref,cord,param,lambda,Qa,Qp,sat,FF, ...
    pos,outvar)
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
%   outvar: output variable, 'x', 'u' or 'x_OL'
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
    
    if strcmp(outvar,'x') || strcmp(outvar,'u')
        uhinf = d2r + 2*lambda*dr+ lambda^2*r;
        u = sat*tanh((Va+Ga+Ma*(-2*lambda*dqa-lambda^2*qa+uhinf))/ ...
            sat); %control signal
    elseif strcmp(outvar,'x_OL')
        u = Va+Ga+Ma*(-2*lambda*dqa-lambda^2*qa); %control signal
    end

end