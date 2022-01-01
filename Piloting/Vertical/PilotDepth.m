function out = PilotDepth(in, memory, parameters)
% Depth piloting
coder.extrinsic('lqrd');
Iz = memory.int_z;
% Iz = EcaF_Saturate(Iz,-parameters.delta_z_sat_m,parameters.delta_z_sat_m);
[BAR0, w0, q0, theta0, A, B, C] = jacob(in);
<<<<<<< Updated upstream
=======

>>>>>>> Stashed changes
dw = in.w_ms - w0;
dq = in.q_rads - q0;
dtheta = in.theta_rad - theta0;
dz = -in.zc_m + in.z_m; 
dz = EcaF_Saturate(dz,-parameters.delta_z_sat_m,parameters.delta_z_sat_m);
<<<<<<< Updated upstream
Q = diag([1e-6,1,1e-6,1e-6,1]);
=======
Q = diag([1,1,1,1,80]);
>>>>>>> Stashed changes
R = 1e8;
Q4 = C.'*Q*C; 
K = zeros(1,5);
K = lqrd(A,B,Q4,R,in.delta_time_s); %wait for teachers solution
<<<<<<< Updated upstream
BARc_rad = (BAR0 + dq*K(1) + dw*K(2) + dz*K(3) + dtheta*K(4) + (Iz - 0)*K(5));

=======
BARc_rad = (BAR0 + dq*K(1) + dw*K(2) + dz*K(3) + dtheta*K(4) + (Iz - 0)*K(5))*0.025;


>>>>>>> Stashed changes
out.BARc = EcaF_Saturate(BARc_rad,parameters.BAR_sat_rad(1),parameters.BAR_sat_rad(2));

end
