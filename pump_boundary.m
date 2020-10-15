function qh_ann = pump_boundary(Qvec,L_eigv,q_pump_initial,area_ann,max_itr,tol)
% Function to calculate the flux at the left boundary
% Newton solver is used for the calculation
% q12np = q1/2 (n+1)
% arguments for newton solver are Function to find solution for, Jacobian
% of the function, initial guess and parameters like maximum number of
% iterations and tolerance.

qh_ann = newton_solver(@(qh_ann)pump_bc(qh_ann,Qvec,L_eigv,q_pump_initial,area_ann),...
        @(qh_ann)pump_bc_jacob(qh_ann,L_eigv),Qvec(:,1),max_itr,tol);

end

