function out = PilotDepth(in, memory, parameters)
% in { u_ms, w_ms, q_rads, theta_rad, z_m , thetac_rad, zc_m , delta_time_s}%
% memory { int_z, current_error, previous_error}%
% params {Kp,Ki,Kd,limit_integral_term,BAR_sat_rad,delta_z_sat_m}%
% Depth piloting
coder.extrinsic('lqrd'); 
%% Compute deltas 
C = eye(5);
Q = eye(5);%diag([1e-6,1,100,10,0.1]);
R = 1;
[BARo, thetao] = smlnkGetEquilibriumPoint(in.u_ms);
wo = 2*tan(thetao); 
qo=0;
[Amat, Bmat] = smlnkGetMatricesAB(in.u_ms, thetao, qo, wo);
 

%% Compute deltas
delta_z = in.zc_m - in.z_m;
delta_z = EcaF_Saturate(delta_z,-parameters.delta_z_sat_m,parameters.delta_z_sat_m);

delta_w = in.w_ms - wo;
delta_q = in.q_rads;
delta_theta = in.theta_rad - thetao;

Qy = C.'*Q*C;

K = zeros(1,5);
K = lqrd(Amat,Bmat,Qy,R,in.delta_time_s);

BARc_rad = BARo + delta_q.*K(1) + delta_w.*K(2) + delta_z.*K(3) + delta_theta.*K(4) + (memory.int_z)*K(5);

%% Output saturation
out.BARc = EcaF_Saturate(BARc_rad,parameters.BAR_sat_rad(1),parameters.BAR_sat_rad(2));
end
