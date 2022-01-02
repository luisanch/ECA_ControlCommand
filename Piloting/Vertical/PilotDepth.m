function out = PilotDepth(in, memory, parameters)
% Depth piloting
coder.extrinsic('lqrd'); 
C = eye(5); 
Q = diag([1e-6,1,100,10,0.0005]);
R = 1e9; 

[BARo,thetao] = smlnkGetEquilibriumPoint(in.u_ms);
[Amat, Bmat] = smlnkGetMatricesAB(in.u_ms, in.theta_rad, in.q_rads , in.w_ms);

qo = 0;
wo = in.u_ms*tan(thetao);

delta_z = in.z_m - in.zc_m;
delta_z = EcaF_Saturate(delta_z,-parameters.delta_z_sat_m,parameters.delta_z_sat_m);
delta_w = in.w_ms - wo;
delta_q = in.q_rads - qo;
delta_theta = in.theta_rad - thetao;

Qy = C'*Q*C;

K = zeros(1,5);
K = lqrd(Amat,Bmat,Qy,R,in.delta_time_s);  
BARc_rad = delta_w.*K(1)+delta_q.*K(2)+delta_z.*K(3)+delta_theta.*K(4) + memory.int_z*K(5); 

%% Output saturation
out.BARc = EcaF_Saturate(BARc_rad,parameters.BAR_sat_rad(1),parameters.BAR_sat_rad(2));

end
