function [BAR_eq, w_eq, q_eq, theta_eq, jacobian_A_double, jacobian_B_double, C] = jacob(in)
% [Cxf, Ksh] = forward_motion();
% [CYuv, CYur, CNuv, CNur] = horizontal();
% [CZuw,CZuq,CMuw,CMuq] = vertical();
CZuw = -3.1430;
CZuq = -1.9597;
CMuw = 0.9709;
CMuq = -0.8872;
u_0 = in.u_ms;
[BAR_eq, theta_eq] = liniarization(u_0);
w_eq = 2*tan(theta_eq);
q_eq = 0;
CM0 = -0.02;
CZ0 = -0.02;
Zb = 0.0;
Zg = 0.03;
CM = -0.17;
m = 1974.26;
CZ = -0.41;
M_inverse = 1e-03*[0.9312 0 0;
    0 0.0514 0;
    0 0 0.2276];



jacobian_A_double = [ 0.0052*CZuw*u_0 (0.0061*q_eq + 0.1015*u_0 + 0.0301 * CZuq*u_0) 0 0 0;
    (0.1333*CMuw*u_0 - 0.0135*q_eq) (0.7758*CMuq*u_0 -0.0135*w_eq) 0 -0.1322*cos(theta_eq) 0
1 0 0 0 0;
0 1 0 0 0;
0 0 1 0 0];

jacobian_B_double = [CZ*u_0^2*0.0052;
    CM*u_0^2*0.1333
    0;
    0;
    0];

C = eye(5);

end

