function l_bc = pump_bc_plant(qh_ann,Qvec_ann,L_eigv,q_pump_initial,area_ann,xi,matA,matB,press_res,in_res,dr,Cl,rho0,p0,rho_ss)
% Function to generate left boundary condition
% The boundary conditions are implicit in q1 and q2 in this case.
% Arguments are the q1/2 (n+1), Q1n and parameters


u = (Cl/2)*qh_ann(1)+(1/2)*qh_ann(2);
v = -(Cl/2)*qh_ann(1)+(1/2)*qh_ann(2);
beta = -rho0*Cl^2 + p0;

rho_ss = qh_ann(1);

% matA(1,1) = ( -1/dr - area_ann/(xi*rho_ss*Cl)); 
% matA(1,2) = (1/dr);
% matC = [-( (area_ann*beta/(xi*Cl*rho_ss)) - (2*area_ann*v/(xi*rho_ss)) + (q_pump_initial)/xi ), zeros(1,length(in_res)+1)]';


% matA(1,1) = ( -1/dr ); 
% matA(1,2) = (1/dr);
% matC = [-( (area_ann/xi)*(qh_ann(2)/qh_ann(1)) + q_pump_initial/xi ) , zeros(1,length(in_res)+1)]';

phi_rw_sub = (qh_ann(1)*Cl^2) + beta;
matC = [phi_rw_sub, zeros(1,length(in_res)+1) ]';

press_res_sub = matA\(matC - matB*press_res);
phi_r_rw = (press_res_sub(2)-press_res_sub(1));
dpdr = phi_r_rw/dr;

q_res_sub = max(xi*dpdr,0);

h1 =  (q_pump_initial+q_res_sub)/area_ann - qh_ann(2)/rho_ss;     % Calculate h1
% h1 = qh_ann(2)/rho_ss - (q_pump_initial+q_res_sub)/area_ann;     % Calculate h1
% h1 = qh_ann(2)/qh_ann(1) - (xi*phi_r_rw + q_pump_initial)/area_ann;
% h1 = (phi_rw_sub-beta)/Cl+2*v - (xi*phi_r_rw/area_ann)*(phi_rw_sub-beta)/Cl^2 - (q_pump_initial/area_ann)*(phi_rw_sub-beta)/Cl^2; %phirrw
h2 = L_eigv(2,:)*((qh_ann-(2*Qvec_ann(:,2)-Qvec_ann(:,3))));     % Calculate h2
% h2 = v - (2*v2-v3);
% h3 = qh_ann(3) - phi_r_rw;% - (phi_rw_sub-beta)/Cl + v - u;% - (qh_ann(1)-rho0)*Cl^2 + p0;

l_bc = [h1;
        h2];                             % Return the boundary condition
    
phi_rw_sub = (qh_ann(1)*Cl^2) + beta;
% putvar(q_res_sub,phi_rw_sub,press_res_sub);

end