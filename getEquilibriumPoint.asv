%function  [bar_eq, theta_eq] = getEquilibriumPoint(u_0)
syms  theta bar u_0

addpath(genpath('.'));
veh_model = loadjson('AUVParameters.json');

m = veh_model.Mecanic.mass_with_on_board_water; 
c_z_uw = veh_model.Hydrodynamic.CZuw; 
c_m_uw = veh_model.Hydrodynamic.CMuw;  
c_z = veh_model.XRearHelms.CZ;
c_m = veh_model.XRearHelms.CM;
rho = veh_model.Environment.Rho;
S_ref = veh_model.Mecanic.surface_reference;
l_ref = veh_model.Mecanic.length_reference;
z_b = 0; 
z_g = 0.03 ;
c_z_0 = -0.02; 
c_m_0 = -0.02;

x1 = (1/2)*l_ref*S_ref*rho*(c_m_0*u_0^2+c_m_uw* theta*u_0^2+bar*c_m*u_0^2)+m*9.81*(z_b-z_g)* theta;
x2 = (1/2)*S_ref*rho*u_0^2*(c_z_0+c_z_uw* theta+bar*c_z);
[bar_eq, theta_eq] = solve([x1==0, x2==0],[bar  theta]); 
% bar_eq = subs(bar_eq);
% theta_eq = subs(theta_eq);
%end
