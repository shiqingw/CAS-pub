function [K_vu,K_vv] = WellboreControlKernels(well_mesh,sig_pp,sig_pm,...
    sig_mp,sig_mm,Cl,q)

well_mesh = well_mesh(2:end-1);

dx = diff(well_mesh);
dx = [dx,dx(end)];
dxi = diff(well_mesh);
dxi = [dxi,dxi(end)];

n = length(well_mesh);

K_uv = nan(n+1,n);
K_vv = nan(n+1,n);

K_vu = nan(n+1,n);
K_uu = nan(n+1,n);

sig_pp = ones(1,n+1)*sig_pp;
sig_pm = ones(1,n+1)*sig_pm;
sig_mp = ones(1,n+1)*sig_mp;
sig_mm = ones(1,n+1)*sig_mm;

%% [K_uv, K_uu ]

for j=1:n
    K_uv(j,j) = sig_pm(j)/(2*Cl);
    K_uv(j+1,j) = sig_pm(j)/(2*Cl); % ghost cells
end

K_uu(1,1) = (1/q)*K_uv(1,1);

for j=2:n
    
    for i=1:j
        
        
        if i <= j-1
            
            K_uv(i,j) = dx(j-1) * ( (sig_pp(j-1)-sig_mm(i))/Cl + 1/dx(j-1)...
                - 1/dxi(i) ) * K_uv(i,j-1) + dx(j-1)/dx(i) * K_uv(i+1,j-1)...
                - (dx(j-1)*sig_pm(i)/Cl )*K_uu(i,j-1);
        end
        
        if i==1
            K_uu(i,j) = (1/q)*K_uv(i,j); % boundary condition
            K_uu(j,j-1) = K_uu(j-1,j-1); % ghost cells
        end

        if i>=2
            K_uu(i,j) = dx(j-1)*( (sig_pp(j-1)-sig_pp(i))/Cl + (1/dx(j-1))...
                - (1/dxi(i)) )*K_uu(i,j-1) + (dx(j-1)/dxi(i))*K_uu(i-1,j-1)...
                - (dx(j-1)*sig_mp(i)/Cl)*K_uv(i,j-1);
            
        end
        
    end
end

%% [K_vu, K_vv]

for j=1:n
    K_vu(j,j) = -sig_mp(j)/(2*Cl);
    K_vu(j+1,j) = -sig_mp(j)/(2*Cl);
end

K_vv(1,1) = (q)*K_vu(1,1);

for j=2:n
    
    for i=1:j
        
        
        if i <= j-1
            
            K_vu(i,j) = dx(j-1) * ( -(sig_mm(j-1)-sig_pp(i))/Cl + 1/dx(j-1)...
                - 1/dxi(i) ) * K_vu(i,j-1) + dx(j-1)/dx(i) * K_vu(i+1,j-1)...
                +(dx(j-1)*sig_mp(i)/Cl )*K_vv(i,j-1);
        end
        
        if i==1
            K_vv(i,j) = q * K_vu(i,j); % boundary condition
            K_vv(j,j-1) = K_vv(j-1,j-1); % ghost cells
        end

        if i>=2
            K_vv(i,j) = dx(j-1)*( -(sig_mm(j-1)-sig_mm(i))/Cl + (1/dx(j-1))...
                - (1/dxi(i)) )*K_vv(i,j-1) + (dx(j-1)/dxi(i))*K_vv(i-1,j-1)...
                + (dx(j-1)*sig_pm(i)/Cl)*K_vu(i,j-1);
            
        end
        
    end
end



end








