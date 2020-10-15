function l_bc = pump_bc(qh_ann,Qvec_ann,L_eigv,q_pump_initial,area_ann)
% Function to generate left boundary condition
% The boundary conditions are implicit in q1 and q2 in this case.
% Arguments are the q1/2 (n+1), Q1n and parameters


h1 = qh_ann(2)/qh_ann(1) - q_pump_initial/area_ann;     % Calculate h1
h2 = L_eigv(2,:)*((qh_ann-(2*Qvec_ann(:,2)-Qvec_ann(:,3))));     % Calculate h2

l_bc = [h1;h2];                             % Return the boundary condition

end