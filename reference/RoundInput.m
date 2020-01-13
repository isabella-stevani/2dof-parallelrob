function [q0,dq0,d2q0] = RoundInput(t)
% Generates 2x1 circular signal with center at (0,0.17), 0.07 radius and
% 2 Hz frequency. Its derivatives are computed as well.
% Inputs:
%   t: time vector
% Outputs:
%   q: circular signal [2x1]
%   dq: first derivative [2x1]
%   d2q: second derivative [2x1]

    C = [0;0.17]; R = 0.07; f = 2; w = 2*pi*f;

    q0 = [C(1) + R*cos(w*t);
          C(2) + R*sin(w*t)];
    
    dq0 = [-w*R*sin(w*t);
           w*R*cos(w*t)];
       
    d2q0 = [-w^2*R*cos(w*t);
           -w^2*R*sin(w*t)];

end