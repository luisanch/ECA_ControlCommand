function [Amat, Bmat] = smlnkGetMatricesAB(u_0, theta, q, w)
Amat = [ -0.1624932128852184532292287855766*u_0, 0.031008020606466744848383313787904*q + 0.010209230070252071084451069569398*u_0, 0, 0, 0;
    0.12682016258403569288959332834007*u_0 - 0.0078338025866952625780154605139582*q, - 0.11588524579777898391886446742436*u_0 - 0.0078338025866952625780154605139582*w, 0, -0.076849603375480526832481280687207*cos(theta), 0;
    1.0, 0, 0, 0, 0;
    0, 1.0, 0, 0, 0;
    0, 0, 1.0, 0, 0];

Bmat = [-0.021197292783548333054112916736472*u_0*abs(u_0);
    -0.022204655633038212801257213664663*u_0*abs(u_0);
    0;
    0;
    0];
end