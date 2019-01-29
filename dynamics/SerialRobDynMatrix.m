function [M,V,G] = SerialRobDynMatrix(q,dq,param)
% Serial mechanism dynamics' matrices.
% Inputs:
%   q: mechanism coordinates
%   dq: mechanism coordinates derivative
%   param: model parameters
% Outputs:
%   M: inertia matrix [2x2]
%   V: centrifugal and coriolis terms vector [2x2]
%   G: gravitational forces vector [2x2]

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

    % RR serial mechanism dynamic model
    D11 = m1*lg1^2 + Jz1 + m2*(l1^2 + lg2^2) + Jz2 + 2*l1*cos(q(2))*m2*lg2;
    D12 = m2*lg2^2 + Jz2 + l1*cos(q(2))*m2*lg2;
    D21 = D12;
    D22 = m2*lg2^2 + Jz2;

    D111 = 0;
    D222 = D111;
    D212 = D111;
    D221 = D111;

    D112 = -l1*sin(q(2))*m2*lg2;
    D121 = D112;
    D122 = D112;

    D211 = l1*sin(q(2))*m2*lg2;
    D1 = g*cos(q(1))*(m1*lg1 + m2*l1) + g*cos(q(1)+q(2))*m2*lg2;
    D2 = g*cos(q(1)+q(2))*m2*lg2;
    
    % Dinamics' matrices
    M = [D11 D12;
        D12 D22];
    V = [D111*dq(1)^2 + D122*dq(2)^2 + D112*dq(1)*dq(2) + D121*dq(1)*dq(2);
        D211*dq(1)^2 + D222*dq(2)^2 + D212*dq(1)*dq(2) + D221*dq(1)*dq(2)];
    G = [D1;D2];

end