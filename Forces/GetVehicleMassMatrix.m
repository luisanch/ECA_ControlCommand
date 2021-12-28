function total_vehicle_mass_matrix = GetVehicleMassMatrix(parameters)

%% Get the parameters and intermediate calculation
m_kg = parameters.mass_with_on_board_water;
zG_m = parameters.center_of_gravity_veh_ref(3);

S = m_kg*[0     zG_m	0;...
          -zG_m 0       0;...
          0     0       0 ];

%% Compute Vehicle "mecaninal mass"
mecanic_matrix= [...
    eye(3)*m_kg    S; ...
    -S             parameters.inertia_matrix]; 

total_vehicle_mass_matrix= mecanic_matrix + parameters.added_water_matrix;
end