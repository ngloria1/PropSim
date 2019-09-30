function mass = CalcMass(results)
psi_to_Pa = 6894.75729; % 1 psi in Pa
    in_to_m = 0.0254; % 1 in in m
    mm_to_m = 1e-3; % 1 mm in m
    lbf_to_N = 4.44822162; % 1 lbf in N
    lbm_to_kg = 0.453592; % 1 lbm in kg
    atm_to_Pa = 101325; % 1 atm in Pa
    L_to_m3 = 1e-3; % 1 L in m^3
    rho_fuel = 870; %kg/m^3
    rho_ox = 1220; %kg/m^3
    mass_engine_cad = 10; %kg
    v_ox_tank_cad = .2;
    v_fuel_tank_cad = .1;
    ox_tank_mass_cad = 8;
    fuel_tank_mass_cad = 5;
    
    v_fuel_tank = results.fuel.V_tank; %m3
    v_ox_tank = results.ox.V_tank; %m3
    fuel_press_init = results.fuel_pressurant.set_pressure; %Pa
    ox_t_init = results.ox.T_tank; %K
    ox_press_init = 630 * psi_to_Pa; %Pa
    
    ox_tank_mass = ox_tank_mass_cad * (v_ox_tank / v_ox_tank_cad);
    fuel_tank_mass = fuel_tank_mass_cad * (v_fuel_tank / v_fuel_tank_cad);
    
    fuel_mass = rho_fuel * v_fuel_tank;
    ox_mass = rho_ox * v_ox_tank;
    
    mass = fuel_mass + ox_mass + fuel_tank_mass + ox_tank_mass + mass_engine_cad;
end