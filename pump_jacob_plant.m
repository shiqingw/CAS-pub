function jacobians = pump_jacob_plant(qh_ann,area_ann,L_eigv,xi,Cl,dr,rho0,p0,q_pump_initial,rho_ss)
% Function to calculate the Jacobian matrix for the left boundary condition
% Left boundary condition is constant volumetric flowrate
% Arguments:
% q12np = q1/2 (n+1)
% param = parameters (structure)


% Get the required parameters
% beta = -rho0*Cl^2+p0;

% phi_r_rw = qh_ann(1);
% phi_rw_sub = qh_ann(3);

% Calculating the terms required for the Jacobian
% J11 = -Cl^2/dr; %-qh_ann(2)/qh_ann(1)^2;%  -(xi/area_ann)*(phi_rw_sub - beta)/Cl^2;
% J11 = -Cl^2/dr + (area_ann/xi)*(qh_ann(2)/qh_ann(1)^2);
J11 = -xi*Cl^2/area_ann;
% J12 = -area_ann/(xi*qh_ann(1));
J12 = -1/rho_ss;
% J13 = -xi/area_ann;%1/Cl - (xi*phi_r_rw)/(area_ann*Cl^2) - q_pump_initial/(area_ann*Cl^2) ;

J21 = L_eigv(2,1);
J22 = L_eigv(2,2);
% J23 = 0;
% 
% J31 = 0;
% J32 = 0;
% J33 = 1;

% J1 = -xi*Cl^2/(area_ann*dr)+ (qh_ann(2)/qh_ann(1)^2);   % J11
% J2 = -(1/qh_ann(1));               % J12
% J3 = L_eigv(2,1);                % J21
% J4 = L_eigv(2,2);                % J22

jacobians = [J11 J12;% J13;
             J21 J22];% J23;
%              J31 J32 J33];  % Return the Jacobian matrix

end
