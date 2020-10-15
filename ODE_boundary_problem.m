function dq = ODE_boundary_problem(x,q,theta,area_ann,dia_bit,dx,g,mu_mud,Cl)
% for source terms
% 2 - drillstring
% 1 - annulus

source_ann = source_terms(q,theta,area_ann,dia_bit,dx,g,mu_mud,1);

dq = [source_ann/Cl^2;
    0];

end