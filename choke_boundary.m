function qh_ann = choke_boundary(Qvec_ann,L_eigv,choke_position,Kc,...
    area_ann,Cl,CTE,max_itr,tol)


qh_ann = newton_solver(@(qh_ann)choke_bc(qh_ann,Qvec_ann,L_eigv,choke_position,Kc,area_ann,Cl,CTE),...
        @(qh_ann)choke_bc_jacob(qh_ann,Kc,choke_position,Cl,CTE,L_eigv,...
        area_ann),Qvec_ann(:,end),max_itr,tol);


    