function [x1] = rk4(x0,u0,h,ffunc)
% Runge-Kutta (4th order).
    K1 = ffunc(x0,u0);
    K2 = ffunc(x0+(h/2)*K1,u0);
    K3 = ffunc(x0+(h/2)*K2,u0);
    K4 = ffunc(x0+h*K3,u0);
    
    x1 = x0+(h/6)*(K1+2*K2+2*K3+K4);
end