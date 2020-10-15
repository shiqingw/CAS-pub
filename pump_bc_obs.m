function l_bc = pump_bc_obs(qh_ann,Qvec_ann,L_eigv,q_pump_initial,area_ann,xi,matA,matB,press_res,in_res,dr,Cl,rho0,p0,BHCP_plant,l1,lr,rho_ss)
% Function to generate left boundary condition
% The boundary conditions are implicit in q1 and q2 in this case.
% Arguments are the q1/2 (n+1), Q1n and parameters


u = (Cl/2)*qh_ann(1)+(1/2)*qh_ann(2);
v = -(Cl/2)*qh_ann(1)+(1/2)*qh_ann(2);
beta = -rho0*Cl^2 + p0;

matA(1,1) = (l1 - 1/dr - 1/(xi*rho_ss*Cl));
matA(1,2) = (1/dr);

phi_rw_sub = qh_ann(1)*Cl^2 + beta;

matC = [-( (beta/(xi*Cl*rho_ss)) - (2*v/(xi*rho_ss)) + (q_pump_initial)/xi - l1*BHCP_plant) , (BHCP_plant-phi_rw_sub).*lr, 0]';

press_res_sub = matA\(matC - matB*press_res);
dpdr = (press_res_sub(2)-press_res_sub(1))/dr;

q_res_sub = max(xi*dpdr,0);

h1 = qh_ann(2)/rho_ss - (q_pump_initial+q_res_sub)/area_ann;     % Calculate h1
h2 = L_eigv(2,:)*((qh_ann-(2*Qvec_ann(:,2)-Qvec_ann(:,3))));     % Calculate h2

l_bc = [h1;
        h2];                             % Return the boundary condition
    


putvar(q_res_sub,phi_rw_sub,press_res_sub);

end