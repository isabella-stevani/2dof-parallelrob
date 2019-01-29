function [q,dq,d2q] = ConstantInput(t)
% Generates 2x1 constant signal with amplitudes 0.1 and 0.2 and its derivatives
% as null vectors.
% Inputs:
%   t: time vector (not used)
% Outputs:
%   q: constant signal [2x1]
%   dq: first derivative [2x1]
%   d2q: second derivative [2x1]

    q = [0.1;0.2];
    dq = zeros(2,1);
    d2q = zeros(2,1);

end