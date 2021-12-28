function out = PilotDepth(in, memory, parameters)
% Depth piloting

%% Compute deltas
delta_z = in.zc_m - in.z_m;
delta_z = EcaF_Saturate(delta_z,-parameters.delta_z_sat_m,parameters.delta_z_sat_m);

delta_w = 0.0 - in.w_ms;

BARc_rad = parameters.Kp*delta_z + parameters.Ki*memory.int_z*delta_z + parameters.Kd*delta_w;

%% Output saturation
out.BARc = EcaF_Saturate(BARc_rad,parameters.BAR_sat_rad(1),parameters.BAR_sat_rad(2));

end
