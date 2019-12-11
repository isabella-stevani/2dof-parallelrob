function [q,dq,u] = ParallelRobDynamics(x0,t,param,uncparam,lambda, ...
    ref,cord,T,sat,FF,loop)
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
%   FF: feedforward enable
%   loop: loop type
% Outputs:
%   q: mechanism coordinates [6x1]
%   dq: mechanism coordinates derivative [6x1]
%   u: control signal [2x1]

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
    
    lambdabar = 10*lambda; %quick convergence for coupling equations
    
    len_t = length(t);
    len_u = 2;
    
    % Feedback linearization matrices previously computed when feedforward
    % is enable
    if FF.on
        Ma = zeros(len_u,len_u,len_t);
        Va = zeros(len_u,len_t);
        Ga = zeros(len_u,len_t);
        
        q = ref(1:6,:);
        dq = ref(7:12,:);
        
        for i = 1:len_t
            [~,~,~,C,dC] = ParallelRobKinMatrix(q(:,i), ...
                dq(:,i),param,Qa,Qp);
            [M,V,G] = ParallelRobDynMatrix(q(:,i), ...
                dq(:,i),param);
            dqa = Qa'*dq(:,i);
            Ma(:,:,i) = C'*M*C;
            Va(:,i) = C'*(M*dC*dqa+V);
            Ga(:,i) = C'*G;
        end
        FF.Ma = Ma;
        FF.Va = Va;
        FF.Ga = Ga;
    end
    
    if strcmp(loop,'CL')
        [~,x] = ode45(@(t,x) ParallelRobDynEDO(t,x,param,uncparam, ...
            lambda,lambdabar,ref,Qa,Qp,cord,T,sat,'x',FF),t,x0);

        x = x';
        q = x(1:6,:);
        dq = x(7:12,:);
        u = zeros(len_u,len_t);

        for i = 1:length(t)
            u(:,i) = ParallelRobDynEDO(t(i),x(:,i),param,uncparam, ...
                lambda,lambdabar,ref,Qa,Qp,cord,T,sat,'u',FF);
        end
    elseif strcmp(loop,'OL')
        [~,x] = ode45(@(t,x) ParallelRobDynEDO(t,x,param,uncparam, ...
            lambda,lambdabar,ref,Qa,Qp,cord,T,sat,'x_OL',FF),t,x0);
        x = x';
        q = x(1:6,:);
        dq = x(7:12,:);
        u = zeros(len_u,len_t);
    end

end
