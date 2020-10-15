function [dx,dt,timesteps,well_mesh,in_mesh,bulk_mud,dia_hyd,area_ann,area_dp,Cl,...
    A_plus,A_minus,L_eigv,sig_pp,sig_pm,sig_mp,sig_mm] = ...
    prepare_well(MD,N_cell,CFL,compress_mud,mu_mud,rho0,run_time,dia_bit,...
    dia_dp_out,dia_dp_in,g)

bulk_mud = 1/compress_mud;
Cl = sqrt(bulk_mud/rho0);
dx = MD/N_cell;
dt = CFL*dx/Cl;
timesteps = floor(run_time/dt);
well_mesh = [0,(0.5*dx):dx:(MD-0.5*dx),MD];
dia_hyd = sqrt(dia_bit^2 - dia_dp_out^2);
area_ann = (pi*dia_hyd^2)/4;
area_dp = (pi*dia_dp_in^2)/4;
in_mesh = 2:length(well_mesh)-1;

A = [0 1; Cl^2 0];
[A_plus,A_minus,L_eigv] = prepare_A(A);

sig_pp = -(1/2)*(32*mu_mud/(rho0*dia_bit^2) + g/Cl);
sig_pm = (1/2)*(-32*mu_mud/(rho0*dia_bit^2) + g/Cl);
sig_mp = -(1/2)*(32*mu_mud/(rho0*dia_bit^2) + g/Cl);
sig_mm = (1/2)*(-32*mu_mud/(rho0*dia_bit^2) + g/Cl);

end