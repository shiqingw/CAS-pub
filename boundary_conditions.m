function res = boundary_conditions(ya,yb,q_pump,area_ann,Cl,choke_position,Kc,rho0,p0,p_atm)
% Function to describe the boundary condition
% Assign values to parameters at both ends
% 0-> L: drillstring & annulus

q1a0 = ya(1);
q2a0 = ya(2);
q1aL = yb(1);
q2aL = yb(2);

res = [q2a0/q1a0 - q_pump/area_ann;
  area_ann^2 * q2aL^2 - (Kc*choke_position)^2 *((2*(q1aL)^2 * Cl^2 - ...
  (2*(rho0 * Cl^2 + (p_atm - p0))) * q1aL))];

end
