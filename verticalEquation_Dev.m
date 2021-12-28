clear variables
clc

import /Forces/.*
%% Parameters
%syms bar  g d a_u_0
%syms rho S_ref  l_ref  W v r m 
% syms m_11_t m_33_t m_55_t m_15 m_35  m_13
% syms det_11 det_33 det_55
% syms c_n c_n_uv c_n_ur
 syms c_x c_x_0 a g  d
% syms c_y c_y_uv c_y_ur
% syms c_z c_z_uw c_z_uq c_z_0
% syms c_m c_m_uw c_m_uq c_m_0
% syms z_b z_g
% syms M_inv
syms theta q bar u_0 w
%% Load path of the workspace
addpath(genpath('.'));
veh_model = loadjson('AUVParameters.json');
M =  GetVehicleMassMatrix(veh_model.Mecanic);

%% 
 %q = 0;
 %w = u_0*tan(0);
 m = veh_model.Mecanic.mass_with_on_board_water;
 W = m * 9.81;
 m_11_t = M(1,1);
 m_33_t = M(3,3);
 m_55_t  = M(5,5);
 m_15 = M(1,5);
 m_35 = M(3,5);
 m_13 = M(1,3);
 c_z_uw = -3.142958771350579;
 c_z_uq = -1.959706445563227;
 c_m_uw = 0.970941769851539;
 c_m_uq = -0.887223477418409;
 c_z = -0.41;
 c_m = -0.17;
 rho = veh_model.Environment.Rho; 
 S_ref = veh_model.Mecanic.surface_reference;
 l_ref = veh_model.Mecanic.length_reference; 
 z_b = 0;
 %u_0 = 2;
 z_g = 0.03 ; 
 c_z_0 = -0.02;
 c_m_0 = -0.02;
 %bar = -0.0896;
 %theta = 0.0053; 


%% Vectors
position_vector = [u_0; w; q];

%% Force matrices
forces_moments = [
    (1/4)*c_x*((a^2) + (bar^2) + (g^2) + (d^2))*u_0*abs(u_0);
    c_z*bar*abs(u_0)*u_0;
    l_ref*c_m*bar*u_0*abs(u_0)];

coriolis = [0, 0 , -m*w;
    0, 0, m*(z_g*q+u_0)
    m*w, -m*(z_g*q+u_0), 0];

damping = [c_x_0, 0, 0;
    0, c_z_uw*u_0, c_z_uq*l_ref*u_0;
    0, c_m_uw*l_ref*u_0, (l_ref^2)*c_m_uq*u_0];

hydrostatic = [0; 0; -(z_g - z_b)*W*sin(theta)];

forces = (1/2)*rho*S_ref*damping*position_vector + coriolis*position_vector + hydrostatic.*position_vector + (1/2)*rho*S_ref*forces_moments;


det_m = m_11_t*m_33_t*m_55_t-m_11_t*(m_35^2)-(m_13^2)*m_55_t+2*m_13*m_15*m_35-(m_15^2)*m_33_t;
det_11 = (m_33_t*m_55_t-m_35^2) / det_m;
det_33 = (m_11_t*m_55_t-(m_15^2)) / det_m;
det_55 = (m_11_t*m_33_t-m_13^2) / det_m;
M_inv = diag([det_11, det_33, det_55]);

velocity_vector = M_inv * forces;

dw = velocity_vector(2);
dq = velocity_vector(3);  

 J = [diff(dw,w) diff(dw,q) 0 diff(dw,theta) 0; diff(dq,w) diff(dq,q) 0 diff(dq,theta) 0; 1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0];
Bsym = [diff(dw,bar);diff(dq,bar);0;0;0];

% Amat = (double(subs(J)));
% Bmat = (double(subs(Bsym)));