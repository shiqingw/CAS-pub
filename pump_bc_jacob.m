function jacobians = pump_bc_jacob(q12np,L_eigv)
% Function to calculate the Jacobian matrix for the left boundary condition
% Left boundary condition is constant volumetric flowrate
% Arguments:
% q12np = q1/2 (n+1)
% param = parameters (structure)


% Get the required parameters

% Calculating the terms required for the Jacobian
term1 = -q12np(2,1)/q12np(1,1)^2;   % J11
term2 = 1/q12np(1,1);               % J12
term3 = L_eigv(2,1);                % J21
term4 = L_eigv(2,2);                % J22

jacobians = [term1 term2; term3 term4];  % Return the Jacobian matrix

end
