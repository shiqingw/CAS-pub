function qh_ann = pump_boundary_plant(Qvec_ann,L_eigv,q_pump_initial,area_ann,max_itr,tol,matA,matB,xi,press_res,in_res,dr,Cl,rho0,p0,rho_ss)
% Function to calculate the flux at the left boundary
% Newton solver is used for the calculation
% q12np = q1/2 (n+1)
% arguments for newton solver are Function to find solution for, Jacobian
% of the function, initial guess and parameters like maximum number of
% iterations and tolerance.
% % 
% uN = (Cl/2)*Qvec_ann(1,:) + (1/2)*Qvec_ann(2,:);
% vN = -(Cl/2)*Qvec_ann(1,:) + (1/2)*Qvec_ann(2,:);
% 
% pr = newton_solver(@(pr)giveF(pr,press_res,in,rmesh,ti,mi,ai,ji,si,ni,vi,theta,BHCP),...
%         @(pr)giveJ(pr,press_res,in,ti,mi,ai,ji,si,ni,vi,theta),press_res,max_itr,tol);

qh_ann = newton_solver(@(qh_ann)pump_bc_plant(qh_ann,Qvec_ann,L_eigv,q_pump_initial,area_ann,xi,matA,matB,press_res,in_res,dr,Cl,rho0,p0,rho_ss),...
        @(qh_ann)pump_jacob_plant(qh_ann,area_ann,L_eigv,xi,Cl,dr,rho0,p0,q_pump_initial,rho_ss),Qvec_ann(:,1),max_itr,tol);

end

