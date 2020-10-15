function [L1,L2,L3,P_plus,P_minus] = WellboreKernels(well_mesh,sig_pp,sig_pm,...
    sig_mp,sig_mm,Cl,q,Kc,area_ann,rho0,p0,p_atm,u,v,q1,q2,choke_position_ss,epsilon)

well_mesh = well_mesh(2:end-1);

dx = diff(well_mesh);
dx = [dx,dx(end)];
dxi = diff(well_mesh);
dxi = [dxi,dxi(end)];

n = length(well_mesh);

P_uv = nan(n,n+1);
P_vv = nan(n,n+1);

P_vu = nan(n,n+1);
P_uu = nan(n,n+1);

sig_pp = ones(1,n+1)*sig_pp;
sig_pm = ones(1,n+1)*sig_pm;
sig_mp = ones(1,n+1)*sig_mp;
sig_mm = ones(1,n+1)*sig_mm;

%% [P_uv, P_vv ]

for j=1:n
    P_uv(j,j) = sig_pm(j)/(2*Cl);
    P_uv(j,j+1) = sig_pm(j)/(2*Cl);
end

P_vv(1,1) = (1/q)*P_uv(1,1);

for j=1:n-1
    
    for i=1:j+1
        
        if i==j
            P_uv(j,i+1) = P_uv(j,i);
        end
        
        if i<j+1
            
            P_uv(j+1,i) = dxi(j)*( (sig_mm(j)-sig_pp(i))/Cl + (1/dxi(j))...
                - (1/dx(i)) )*P_uv(j,i) + (dxi(j)/dx(i))*P_uv(j,i+1)...
                - (dxi(j)*sig_pm(i)/Cl )*P_vv(j,i);
            
        end
        
        if i==1
            P_vv(j+1,i) = (1/q)*P_uv(j+1,i);
        end
        
        if i>1
            
            if i==j+1
                P_vv(j,i) = P_vv(j,i-1);
            end
            
            P_vv(j+1,i) = dxi(j)*( (sig_mm(j)-sig_mm(i))/Cl + (1/dxi(j))...
                - (1/dx(i)) )*P_vv(j,i) + (dxi(j)/dx(i))*P_vv(j,i-1)...
                - (dxi(j)*sig_mp(i)/Cl)*P_uv(j,i);
            
        end
        
    end
end

%% [P_vu, P_uu]

for j=1:n
    P_vu(j,j) = -sig_mp(j)/(2*Cl);
    P_vu(j,j+1) = -sig_mp(j)/(2*Cl);
end

P_uu(1,1) = (q)*P_vu(1,1);

for j=1:n-1
    
    for i=1:j+1
        
        if i==j
            P_vu(j,i+1) = P_vu(j,i);
        end
        
        if i<j+1
            
            P_vu(j+1,i) = dxi(j)*( (sig_mm(i)-sig_pp(j))/Cl + (1/dxi(j))...
                - (1/dx(i)) )*P_vu(j,i) + (dxi(j)/dx(i))*P_vu(j,i+1)...
                + (dxi(j)*sig_mp(i)/Cl)*P_uu(j,i);
            
        end
        
        if i==1
            P_uu(j+1,i) = (q)*P_vu(j+1,i);
        end
        
        if i>1
            
            if i==j+1
                P_uu(j,i) = P_uu(j,i-1);
            end
            
            P_uu(j+1,i) = dxi(j)*( (sig_pp(i)-sig_pp(j))/Cl + (1/dxi(j))...
                - (1/dx(i)) )*P_uu(j,i) + (dxi(j)/dx(i))*P_uu(j,i-1)...
                + (dxi(j)*sig_pm(i)/Cl)*P_vu(j,i);
            
        end
        
    end
end

alpha = Cl;
beta = -Cl;

a_t = Kc*choke_position_ss;
b = (-rho0*Cl^2 + p0 - p_atm)*2;

h1 = -(q2/(q1^2))*area_ann + (a_t/sqrt( 2*Cl^2 + b/q1 ))* (b/(q1^2));
h2 = area_ann/q1;

H_bar = (h1+Cl*h2)/(h1-Cl*h2);

Theta = (H_bar - epsilon*alpha)/(1 + epsilon*beta);

P_plus = Cl*( Theta*P_uv(end,1:end-1) - P_uu(end,end-1) )/( alpha+ Theta*beta );
P_minus = Cl*( Theta*P_vv(end,1:end-1) - P_vu(end,end-1) )/( alpha+ Theta*beta );

L1 = (P_plus - P_minus)./Cl;
L2 = P_plus + P_minus;

L3 = (Theta*( h2 - h1/Cl ) + (h2 + h1/Cl))/ ( Cl*(Theta-1) );

% L3 = epsilon*dhdv/ (1+ epsilon*beta + epsilon*Cl);

%%
% aa = Kc*sqrt(2/Cl);
% bb = -rho0*Cl^2 + p0 - p_atm;
% 
% dhdu = area_ann - (aa*choke_position_ss/2)*( 2*(u-v)*Cl + bb )...
%     /sqrt( (u-v)^2 *Cl + (u-v)*bb );
% dhdv = area_ann + (aa*choke_position_ss/2)*( 2*(u-v)*Cl + bb )...
%     /sqrt( (u-v)^2 *Cl + (u-v)*bb );
% 
% ro = -dhdu/dhdv;
% 
% alpha = Cl;
% beta = -Cl;
% 
% P_plus = (ro-epsilon*alpha)*P_uv(end,1:end-1) - ...
%     (1+epsilon*beta)*P_uu(end,1:end-1)*(Cl/(alpha*(1+epsilon*beta) ...
%     + beta*(ro-epsilon*alpha)));
% 
% P_minus = (ro-epsilon*alpha)*P_vv(end,1:end-1) - ...
%     (1+epsilon*beta)*P_vu(end,1:end-1)*(Cl/(alpha*(1+epsilon*beta)...
%     + beta*(ro-epsilon*alpha)));
% % 
% % P_plus = (k_bar*P_uv(end,1:end-1) - P_uu(end,1:end-1))/(1-k_bar);
% % P_minus = (k_bar*P_vv(end,1:end-1) - P_vu(end,1:end-1))/(1-k_bar);
% 
% L1 = (P_plus - P_minus)./Cl;
% L2 = P_plus + P_minus;
% 
% L3 = epsilon*dhdv/ (1+ epsilon*beta + epsilon*Cl);
% % L3 = ((ro-epsilon*beta)*dhdv + (1+epsilon*beta)*dhdu)/(Cl*(ro-1));


end








