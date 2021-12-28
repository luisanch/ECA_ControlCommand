load('../Data/PitchSteps_NoiseOn_Pitch_W_Q.mat');

%data values

Fz = [forces_values.HydrodynamicForces.Fz_N.Data(:)];
My = [forces_values.HydrodynamicForces.My_Nm.Data(:)];
u_ms = vehicle_state.nu.rAUV_WaterSpeed.Uwater_ms.Data;
w_ms = vehicle_state.nu.rAUV_WaterSpeed.Wwater_ms.Data;
q_rads = vehicle_state.nu.AngularSpeed.Q_rads.Data;

rho = veh_model.Environment.Rho;
L_ref = veh_model.Mecanic.length_reference;
S_ref = veh_model.Mecanic.surface_reference;

coef = 0.5*rho*S_ref;

Vertical_force_parameters = (coef.*[ u_ms.*w_ms L_ref.*u_ms.*q_rads])\Fz;
CZuw = Vertical_force_parameters(1);
CZuq = Vertical_force_parameters(2);

Vertical_moment_parameters = (coef.*[ L_ref.*u_ms.*w_ms (L_ref.^2).*u_ms.*q_rads])\My;
CMuw = Vertical_moment_parameters(1);
CMuq = Vertical_moment_parameters(2);