syms uo w q theta rho Sref Lref A BAR Cx G D Cz Cm m zg q Cxo Czuw Cmuw Czuq Cmuq zb W CZo CMo

M1 = [(1/4)*Cx*(A^2+BAR^2+G^2+D^2)*uo*abs(uo);Cz*BAR*uo*abs(uo);Lref*Cm*BAR*uo*abs(uo)];
M2 = [0 0 -m*w; 0 0 m*(zg*q+uo); m*w -m*(zg*q+uo) 0];
M3 = [Cxo 0 0; 0 Czuw*uo Czuq*Lref*uo; 0 Cmuw*Lref*uo Lref*Lref*Cmuq*uo];
M4 = [0; 0; -(zg-zb)*W*sin(theta)];

V = [uo; w; q];

K = (1/2)*rho*Sref*M1+M2*V+(1/2)*rho*Sref*M3*V+M4.*V;

K(2) = K(2)+(rho*Sref*CZo*uo^2)/2;
K(3) = K(3)+(rho*Sref*CMo*uo^2)/2;

syms m11t m33t m55t m13 m15 m35 

DetM = m11t*m33t*m55t-m11t*m35^2-m13^2*m55t+2*m13*m15*m35-m15^2*m33t;
Det11 = (m33t*m55t-m35^2)/DetM;
Det33 = (m11t*m55t-m15^2)/DetM;
Det55 = (m11t*m33t-m13^2)/DetM;

Mm = [Det11 0 0; 0 Det33 0; 0 0 Det55];



R=Mm*K;

dw=R(2);
dq=R(3);

J = [diff(dw,w) diff(dw,q) 0 diff(dw,theta) 0; diff(dq,w) diff(dq,q) 0 diff(dq,theta) 0; 1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0];