function [q1_ann_ss,q2_ann_ss] = CalculateSteadyState(q_pump_ss,area_ann,...
    Cl,Kc,choke_position_ss,well_mesh,in_mesh,dx,theta,...
    dia_bit,g,mu_mud,rho0,p0,p_atm,A_plus,A_minus,L_eigv,dt,time_initial,max_itr,tol)

Qvec_int = [rho0 + (p_atm - p0)/Cl^2;
    (rho0 + (p_atm - p0)/Cl^2)*q_pump_ss/area_ann];

initial_solution = bvpinit(well_mesh,Qvec_int);

solution = bvp4c(@(x,Qvec)ODE_boundary_problem(x,Qvec,theta,area_ann,dia_bit,dx,g,mu_mud,Cl),...
        @(ya,yb)boundary_conditions(ya,yb,q_pump_ss,area_ann,Cl,choke_position_ss,Kc,rho0,p0,p_atm),...
        initial_solution);

y = deval(solution,well_mesh);

Qvec = y(1:2,:);

CTE = 2*(rho0 * Cl^2 + (p_atm - p0));


for t=1:floor(time_initial/dt)
        
%     [Qvec(:,1)] = pump_boundary_plant(Qvec,L_eigv,q_pump_ss,area_ann,max_itr,tol,matA,matB,xi,press_res,in_res,dr,Cl,rho0,p0,rho_ss);
    Qvec(:,1) = pump_boundary(Qvec,L_eigv,q_pump_ss,area_ann,max_itr,tol); 
    Qvec(:,end) = choke_boundary(Qvec,L_eigv,choke_position_ss,Kc,area_ann,Cl,CTE,max_itr,tol);
    
    [F_plus_ann, F_minus_ann] = flux_vec(Qvec, A_plus , A_minus);
    
    Qvec(:,in_mesh) = Qvec(:,in_mesh) - (dt/dx)*(F_plus_ann - F_minus_ann)+...
        dt*[zeros(1,length(in_mesh));source_terms(Qvec(:,in_mesh),theta,area_ann,dia_bit,dx,g,mu_mud,1)];     
    
end

q1_ann_ss = Qvec(1,:);
q2_ann_ss = Qvec(2,:);








