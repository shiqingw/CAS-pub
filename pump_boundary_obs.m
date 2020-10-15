function qh_ann = pump_boundary_obs(Qvec_ann,L_eigv,q_pump_initial,area_ann,max_itr,tol,matA,matB,xi,press_res,in_res,dr,Cl,rho0,p0,BHCP_plant,l1,lr,rho_ss)
% Function to calculate the flux at the left boundary
% Newton solver is used for the calculation
% q12np = q1/2 (n+1)
% arguments for newton solver are Function to find solution for, Jacobian
% of the function, initial guess and parameters like maximum number of
% iterations and tolerance.

qh_ann = newton_solver(@(qh_ann)pump_bc_obs(qh_ann,Qvec_ann,L_eigv,q_pump_initial,area_ann,xi,matA,matB,press_res,in_res,dr,Cl,rho0,p0,BHCP_plant,l1,lr,rho_ss),...
        @(qh_ann)pump_jacob_obs(qh_ann,area_ann,L_eigv,xi,Cl,dr,rho_ss,l1),Qvec_ann(:,1),max_itr,tol);

end

