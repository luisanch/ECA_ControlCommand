function out = PilotDepth(in, memory, parameters)
% Depth piloting
coder.extrinsic('lqrd');
Iz = memory.int_z;
% Iz = EcaF_Saturate(Iz,-parameters.delta_z_sat_m,parameters.delta_z_sat_m);
[BAR0, w0, q0, theta0, A, B, C] = jacob(in);
dw = in.w_ms - w0;
dq = in.q_rads - q0;
dtheta = in.theta_rad - theta0;
dz = -in.zc_m + in.z_m; 
dz = EcaF_Saturate(dz,-parameters.delta_z_sat_m,parameters.delta_z_sat_m);
Q = diag([1e-6,1,1e-6,1e-6,1]);
R = 1e8;
Q4 = C.'*Q*C; 
K = zeros(1,5);
K = lqrd(A,B,Q4,R,in.delta_time_s); %wait for teachers solution
BARc_rad = (BAR0 + dq*K(1) + dw*K(2) + dz*K(3) + dtheta*K(4) + (Iz - 0)*K(5));

% %% Compute deltas
% delta_z = in.zc_m - in.z_m;
% delta_z = EcaF_Saturate(delta_z,-parameters.delta_z_sat_m,parameters.delta_z_sat_m);

% delta_w = 0.0 - in.w_ms;
% 
% BARc_rad = parameters.Kp*delta_z + parameters.Ki*memory.int_z*delta_z + parameters.Kd*delta_w;
% 
% %% Output saturation
% out.BARc = EcaF_Saturate(BARc_rad,parameters.BAR_sat_rad(1),parameters.BAR_sat_rad(2));
out.BARc = EcaF_Saturate(BARc_rad,parameters.BAR_sat_rad(1),parameters.BAR_sat_rad(2));

end
