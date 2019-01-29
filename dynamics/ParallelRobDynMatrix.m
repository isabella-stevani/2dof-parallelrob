function [M,V,G] = ParallelRobDynMatrix(q,dq,param)
% Parallel mechanism dynamics' matrices.
% Inputs:
%   q: mechanism coordinates
%   dq: mechanism coordinates derivative
%   param: model parameters
% Outputs:
%   M: inertia matrix [6x6]
%   V: centrifugal and coriolis terms vector [6x1]
%   G: gravitational forces vector [6x1]

    % Dynamics' matrices for each serial subsystem
    q1 = q(3:4); q2 = q(5:6);
    dq1 = dq(3:4); dq2 = dq(5:6);
    [M1,V1,G1] = SerialRobDynMatrix(q1,dq1,param);
    [M2,V2,G2] = SerialRobDynMatrix(q2,dq2,param);
    
    % Concatenating serial subsystems's matrices
    M = blkdiag(zeros(size(M1)),M1,M2);
    V = [zeros(size(V1));V1;V2];
    G = [zeros(size(G1));G1;G2];

end