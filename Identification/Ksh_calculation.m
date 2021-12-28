load('../Data/SpeedSteps_Noise_U.mat');

%data values

u_ms = vehicle_state.nu.rAUV_WaterSpeed.Uwater_ms.Data;
Fx = [forces_values.HydrodynamicForces.Fx_N.Data(:)];

rho = veh_model.Environment.Rho;
viscosity_nu = veh_model.Environment.nu;
S_ref = veh_model.Mecanic.surface_reference;
L_ref = veh_model.Mecanic.length_reference;

%Calculate Ksh

f_formula = Fx./(0.5*rho*S_ref.*u_ms.*abs(u_ms));

Re = u_ms*L_ref/viscosity_nu;
Cxf = 0.075./(log10(Re)-2).^2;

Ksh = Cxf\f_formula;