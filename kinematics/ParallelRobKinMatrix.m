function [qbar,A,dA,C,dC] = ParallelRobKinMatrix(q,dq,param,Qa,Qp)
% Parallel mechanism dynamics' matrices.
% Inputs:
%   q: mechanism coordinates
%   dq: mechanism coordinates derivative
%   param: model parameters
%   Qa: active coordinates matrix
%   Qp: passive coordinates matrix
% Outputs:
%   qbar: mechanical coupling vector [6x1]
%   A:  coupling matrix A [6x6]
%   dA: A matrix derivative [6x6]
%   C: coupling matrix C [6x2]
%   dC: C matrix derivative [6x2]

    % Model parameters and coordinates
    g = param.g;
    m1 = param.m1;
    m2 = param.m2;
    l0 = param.l0;
    l1 = param.l1;
    l2 = param.l2;
    lg1 = param.lg1;
    lg2 = param.lg2;
    Jz1 = param.Jz1;
    Jz2 = param.Jz2;
    
    q0 = q(1:2); q1 = q(3:4); q2 = q(5:6);
    dq0 = dq(1:2); dq1 = dq(3:4); dq2 = dq(5:6);
    
    qbar = [q0(1)-l0-l1*cos(q1(1))-l2*cos(q1(1)+q1(2));
            q0(2)-l1*sin(q1(1))-l2*sin(q1(1)+q1(2));
            q0(1)+l0+l1*cos(q2(1))+l2*cos(q2(1)+q2(2));
            q0(2)-l1*sin(q2(1))-l2*sin(q2(1)+q2(2))]; %coupling vector
        
    A0 = [1 0;
          0 1;
          1 0;
          0 1];
      
    A1 = [l1*sin(q1(1))+l2*sin(q1(1)+q1(2)) l2*sin(q1(1)+q1(2));
          -l1*cos(q1(1))-l2*cos(q1(1)+q1(2)) -l2*cos(q1(1)+q1(2));
          0 0;
          0 0];
    A2 = [0 0;
          0 0;
          -l1*sin(q2(1))-l2*sin(q2(1)+q2(2)) -l2*sin(q2(1)+q2(2));
          -l1*cos(q2(1))-l2*cos(q2(1)+q2(2)) -l2*cos(q2(1)+q2(2))];
    
    A = [A0 A1 A2]; %dqbar/dq
    
    C = Qa - Qp*((A*Qp)\(A*Qa)); %vector q from active coordinates vector qa
    
    dA0 = zeros(4,2);
    
    dA1 = [l1*cos(q1(1))*dq1(1)+l2*cos(q1(1)+q1(2))*(dq1(1)+dq1(2)) ...
        l2*cos(q1(1)+q1(2))*(dq1(1)+dq1(2));
           l1*sin(q1(1))*dq1(1)+l2*sin(q1(1)+q1(2))*(dq1(1)+dq1(2)) ...
        l2*sin(q1(1)+q1(2))*(dq1(1)+dq1(2));
           0 0;
           0 0];
    
    dA2 = [0 0;
           0 0;
           -l1*cos(q2(1))*dq2(1)-l2*cos(q2(1)+q2(2))*(dq2(1)+dq2(2)) ...
        -l2*cos(q2(1)+q2(2))*(dq2(1)+dq2(2));
           l1*sin(q2(1))*dq2(1)+l2*sin(q2(1)+q2(2))*(dq2(1)+dq2(2)) ...
        l2*sin(q2(1)+q2(2))*(dq2(1)+dq2(2))];
       
    dA = [dA0 dA1 dA2]; %dA/dq
    
    dC = -Qp*((A*Qp)\(dA*C)); %dC/dq

end
