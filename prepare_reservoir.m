function [r_mesh,dr,xi,compress_total,bulk_oil,eps]  = ...
    prepare_reservoir(r_w,r_e,r_cell,compress_oil,compress_res,mu_oil,...
    p0,rho0,phi,perm,h_res,dt,theta_scheme)

r_mesh = [linspace(r_w,r_e,r_cell) r_e+0.1]; % 1 ghost cell at end

compress_total = compress_oil + compress_res;
bulk_oil = 1/compress_oil;
area_flow = 2*pi*r_w*h_res;
xi = perm*area_flow/mu_oil;
eps = perm/(mu_oil*phi*compress_total);


in_res = 2:(length(r_mesh)-1);     % inner domain

rip1 = r_mesh(in_res+1);         % r_(i+1)
ri = r_mesh(in_res);             % r_i
rim1 = r_mesh(in_res-1);         % r_(i-1)

rphalf = (ri + rip1)/2;     % r_(i+1/2)
rmhalf = (ri + rim1)/2;     % r_(i-1/2)
% 
% alpha = rip1 - ri;          % r_(i+1) - r_i
% beta = ri-rim1;             % r_i - r_(i-1)
% gamma = rphalf - rmhalf;    % r_(i+1/2) - r_(i-1/2)
% 
% rphi = ((perm*rho0)/(bulk_oil*mu_oil)).*rphalf; 
% rmhi = ((perm*rho0)/(bulk_oil*mu_oil)).*rmhalf;
% ti = (rho0*phi*compress_total/bulk_oil).*ri;
% mi = (alpha.*beta.*gamma)./dt;
% ai = (bulk_oil-p0)*ti;
% ji = beta.*(bulk_oil-p0).*rphi;
% si = (beta.*rphi)./2; 
% ni = alpha.*(bulk_oil-p0).*rmhi;
% vi = alpha.*rmhi./2;
% 
% putvar(ti,mi,ai,ji,si,ni,vi)

dr = (r_mesh(2)-r_mesh(1));

alphai = ri./(eps*dt);
betai = -rphalf./((rphalf-rmhalf).*(rip1-ri));
gammai = rmhalf./((rphalf-rmhalf).*(ri-rim1));

A1 = -theta_scheme.*gammai;
A2 = alphai + theta_scheme.*gammai - theta_scheme.*betai;
A3 = theta_scheme.*betai;

B1 = -(1-theta_scheme).*gammai;
B2 = -alphai + (1-theta_scheme).*gammai - (1-theta_scheme).*betai;
B3 = (1-theta_scheme).*betai;

matA = diag([A1,-1],-1) + diag([ 1,A2,1]) + diag([ 0,A3],1);
matB = diag([B1,0],-1) + diag([0,B2,0]) + diag([0,B3],1);

putvar(matA,matB,in_res,dr); 

end