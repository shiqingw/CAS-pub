function [dq1_ann,dq2_ann,Qvec] = initialization(q_pump_initial,area_ann,...
    Cl,Kc,choke_position_initial,well_mesh,in_mesh,dx,theta,dia_bit,g,mu_mud,rho0,p0,...
    p_atm,A_plus,A_minus,L_eigv,dt,time_initial,q1_ann_ss,q2_ann_ss,max_itr,tol)

Qvec_int = [rho0 + (p_atm - p0)/Cl^2; 
    (rho0 + (p_atm - p0)/Cl^2)*q_pump_initial/area_ann];

initial_solution = bvpinit(well_mesh,Qvec_int);

solution = bvp4c(@(x,Qvec)ODE_boundary_problem(x,Qvec,theta,area_ann,dia_bit,dx,g,mu_mud,Cl),...
        @(ya,yb)boundary_conditions(ya,yb,q_pump_initial,area_ann,Cl,choke_position_initial,Kc,rho0,p0,p_atm),...
        initial_solution);

y = deval(solution,well_mesh);

Qvec = y(1:2,:);


CTE = 2*(rho0 * Cl^2 + (p_atm - p0));
% 
% u = (Cl/2)*Qvec(1,:)+(1/2)*Qvec(2,:);
% v = -(Cl/2)*Qvec(1,:)+(1/2)*Qvec(2,:);
% 
% matA = [zeros(1,length(in_res)+3);zeros(length(in_res)+2,1),matA];
% 
% matA(1,1) = 1;
% matA(1,2) = -1/Cl;
% 
% matA(2,2) = ( -1/dr - area_ann/(xi*rho_ss*Cl)); 
% matA(2,3) = (1/dr);
% 
% 
% matB = [zeros(1,length(in_res)+3);zeros(length(in_res)+2,1),matB];
% matC = zeros(length(in_res)+3,1);
% 
% beta = -rho0*Cl^2 + p0;
% 
% trysol = zeros(length(in_res)+3,floor(time_initial/dt)+1);
% trysol(:,1) = [u(1);press_res];
% 
% for t=1:floor(time_initial/dt)
% 
%     matC(1,1) = -beta/Cl+v(1);
%     matC(2,1) = 2*area_ann*v(1)/(xi*rho_ss)-q_pump_initial/xi -...
%         area_ann*beta/(rho_ss*xi*Cl);
%     
%     trysol(:,t+1) = matA\(matC-matB*trysol(:,1));
%     
%     dpdr = (trysol(3,t+1)-trysol(2,t+1))/dr;
%     q_res = xi*dpdr
%     
% %     v_t = v(1);
%     v_t = 2*v(2)-v(3);
% %     phirw = trysol(1,t+1);
%     u_t = trysol(1,t);%v_t + (phirw-beta)/Cl;
%     
%     Qvec(1,1) = (u_t-v_t)/Cl;
%     Qvec(2,1) = (u_t+v_t);
% 
% %     Qvec(:,1) = pump_boundary_plant(Qvec,L_eigv,q_pump_initial,...
% %      area_ann,max_itr,tol,matA,matB,xi,press_res,in_res,dr,Cl,rho0,p0,rho_ss);
% %     Qvec(:,1) = pump_boundary(Qvec,L_eigv,q_pump_initial,area_ann,max_itr,tol); 
%     Qvec(:,end) = choke_boundary(Qvec,L_eigv,choke_position_initial,...
%         Kc,area_ann,Cl,CTE,max_itr,tol);
%     
%     [F_plus_ann, F_minus_ann] = flux_vec(Qvec, A_plus , A_minus);
%     
%     Qvec(:,in_mesh) = Qvec(:,in_mesh) - (dt/dx)*(F_plus_ann - F_minus_ann)+...
%         dt*[zeros(1,length(in_mesh));source_terms(Qvec(:,in_mesh),...
%         theta,area_ann,dia_bit,dx,g,mu_mud,1)];    
%     
%     v = -(Cl/2)*Qvec(1,:)+(1/2)*Qvec(2,:);
% 
% end
% press_sub = trysol(2:end,end);
% putvar(press_sub);

for t=1:floor(time_initial/dt)
        
%     [Qvec(:,1)] = pump_boundary_plant(Qvec,L_eigv,q_pump_ss,area_ann,max_itr,tol,matA,matB,xi,press_res,in_res,dr,Cl,rho0,p0,rho_ss);
    Qvec(:,1) = pump_boundary(Qvec,L_eigv,q_pump_initial,area_ann,max_itr,tol); 
    Qvec(:,end) = choke_boundary(Qvec,L_eigv,choke_position_initial,Kc,area_ann,Cl,CTE,max_itr,tol);
    
    [F_plus_ann, F_minus_ann] = flux_vec(Qvec, A_plus , A_minus);
    
    Qvec(:,in_mesh) = Qvec(:,in_mesh) - (dt/dx)*(F_plus_ann - F_minus_ann)+...
        dt*[zeros(1,length(in_mesh));source_terms(Qvec(:,in_mesh),theta,area_ann,dia_bit,dx,g,mu_mud,1)];     
    
end

dq1_ann = Qvec(1,:) - q1_ann_ss;
dq2_ann = Qvec(2,:) - q2_ann_ss;





