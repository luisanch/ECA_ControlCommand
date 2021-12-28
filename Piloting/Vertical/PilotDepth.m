function out = PilotDepth(in, memory, parameters)
% in { u_ms, w_ms, q_rads, theta_rad, z_m , thetac_rad, zc_m , delta_time_s}%
% memory { int_z, current_error, previous_error}%
% params {Kp,Ki,Kd,limit_integral_term,BAR_sat_rad,delta_z_sat_m}%
% Depth piloting
coder.extrinsic('lqrd');

%% Compute deltas
Amat = [
 -0.1624*in.u_ms, 0.0310*in.q_rads + 0.0102*in.u_ms, 0, 0, 0;
 0.1268* in.u_ms - 0.0078*in.q_rads, - 0.5794* in.u_ms - 0.0078*in.w_ms - 0.0768*sin(in.theta_rad), 0, -0.0768*in.q_rads*cos(in.theta_rad), 0;
 1, 0, 0, 0, 0;
 0, 1, 0, 0, 0;
 0, 0, 1, 0, 0
 ];
Bmat = [-0.0211* in.u_ms*abs(in.u_ms); -0.0222*in.u_ms*abs(in.u_ms); 0; 0; 0];

C = eye(5);
Q =diag([1,1,1,1,1]);
R = 1;

BARo = -0.0896;
wo = 2*tan(0.0053);
qo = 0;
thetao = 0.0053; 

%% Compute deltas
delta_z = in.zc_m - in.z_m;
delta_z = EcaF_Saturate(delta_z,-parameters.delta_z_sat_m,parameters.delta_z_sat_m);

delta_w = in.w_ms - wo;
delta_q = in.q_rads - qo;
delta_theta = in.theta_rad - thetao;

Qy = C'*Q*C;

K = zeros(1,5);
K = lqrd(Amat,Bmat,Qy,R,in.delta_time_s);

BARc_rad = BARo + delta_w.*K(1)+delta_q.*K(2)+delta_z.*K(3)+delta_theta.*K(4)+memory.int_z.*K(5); 

%% Output saturation
out.BARc = EcaF_Saturate(BARc_rad,parameters.BAR_sat_rad(1),parameters.BAR_sat_rad(2));

end
