function jacobians = pump_jacob_obs(qh_ann,area_ann,L_eigv,xi,Cl,dr,rho_ss,l1)
% Function to calculate the Jacobian matrix for the left boundary condition
% Left boundary condition is constant volumetric flowrate
% Arguments:
% q12np = q1/2 (n+1)
% param = parameters (structure)


% Get the required parameters

% Calculating the terms required for the Jacobian
J1 = l1*Cl^2 - Cl^2/dr  ;%-xi*Cl^2/(area_ann*dr)+ (qh_ann(2)/qh_ann(1)^2);   % J11
J2 = -1/(xi*rho_ss) ;%-(1/qh_ann(1));               % J12
J3 = L_eigv(2,1);                % J21
J4 = L_eigv(2,2);                % J22

jacobians = [J1 J2;
             J3 J4];  % Return the Jacobian matrix

end
