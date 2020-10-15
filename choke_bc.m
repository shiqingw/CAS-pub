function q_rbc = choke_bc(qh_ann,Qvec_ann,L_eigv,choke_position,Kc,area_ann,Cl,CTE)

h1 = L_eigv(1,:)*(qh_ann-(2*Qvec_ann(:,end-1)-Qvec_ann(:,end-2)));    % h1 for the solver
h2 = area_ann * qh_ann(2) -...          
     sign(2*(qh_ann(1))^2 * Cl^2 - CTE * qh_ann(1))*Kc*choke_position *...
     sqrt(abs(2*(qh_ann(1))^2 * Cl^2 - CTE * qh_ann(1))); % h2 for the solver
 
q_rbc = [h1;h2];                    % return the vector for the solver
