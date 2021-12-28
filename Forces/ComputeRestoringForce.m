function restoring_forces = ComputeRestoringForce(EulerAngles,gravity_center_m,buoyancy_center_m,mass_kg)
% Compute coriolis forces
%
% :usage: coriolis_forces = ComputeCoriolisForces(speed_ms,angular_speed_rads ,...gravity_center_m, mass_kg,inertia_matrix)
%
% :param EulerAngles: [rall,pitch,yaw] in rad
% :param gravity_center_m: 
%   Position of the center of gravity in relation to the center 
%   of expression of the forces
% :param buoyancy_center_m: 
%   Position of the center of buoyancy in relation to the center 
%   of expression of the forces
% :param mass_kg: vehicule mass in Kg
%
% :returns:
%   * restoring_forces - vector of restoring forces [Fx,Fy,Fz,Mx,My,Mz]

% Variables recovery
roll_rad = EulerAngles.Phi_rad;
pitch_rad = EulerAngles.Theta_rad;
W = mass_kg*9.81;

sth  = sin(pitch_rad); cth  = cos(pitch_rad);sphi = sin(roll_rad);

% ATTENTION aux espaces avant et après les parenthèses !
restoring_forces = -[0;...
    0;...
    0;...
    (gravity_center_m(3)-buoyancy_center_m(3))*W*cth*sphi;...
    (gravity_center_m(3)-buoyancy_center_m(3))*W*sth;...
   0];
