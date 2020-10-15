function q_rbc = choke_bc_obs1(qh_ann,Qvec_ann,L_eigv,choke_position,Kc,area_ann,...
    Cl,CTE,pchoke_plant,L3,rho0,p0, C1, C2, K_vu, K_vv, choke_position_ss)

pchoke_obs = (qh_ann(1)-rho0)*Cl^2+p0;

h1 = L_eigv(1,:)*(qh_ann-(2*Qvec_ann(:,end-1)-Qvec_ann(:,end-2)));    % h1 for the solver
h2 = area_ann * qh_ann(2) -...          
     sign(2*(qh_ann(1))^2 * Cl^2 - CTE * qh_ann(1))*Kc*choke_position *...
     sqrt(abs(2*(qh_ann(1))^2 * Cl^2 - CTE * qh_ann(1))) - L3*(pchoke_plant-pchoke_obs); % h2 for the solver
 
q_rbc = [h1;h2];                     % return the vector for the solver

end