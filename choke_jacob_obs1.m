function choke_jacobians = choke_jacob_obs1(qh_ann,Kc,choke_position,...
    Cl,CTE,L_eigv,area_ann,L3,C1, C2, K_vu, K_vv, choke_position_ss)

% Calculate the terms for the Jacobian matrix
J11 = L_eigv(1,1);                % J11
J12 = L_eigv(1,2);                % J12
J21 = - Kc*choke_position * (4 * qh_ann(1) * Cl^2 - CTE)/ ...
   (2*sqrt(abs(2*(qh_ann(1))^2 * Cl^2 - CTE * qh_ann(1)))) + L3*Cl^2;   %J21
J22 = area_ann;                                                   % J22

choke_jacobians = [J11 J12; 
                   J21 J22]; 

end