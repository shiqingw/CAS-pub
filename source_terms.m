function S = source_terms(Qvec,theta,area,dia,dx,g,mu_mud,n)
%Function to calculate the source terms
% n = 1 : annulus
% n = 2 : drill string


theta = (theta/180)*pi;         % Inclination in radians
vel = Qvec(2,:)./Qvec(1,:);             % velocity, m/s

p = length(Qvec(1,:));
F = - Kelda((vel*area)', dia*ones(p,1), area*ones(p,1),...
    dx*ones(p,1), (Qvec(1,:))', mu_mud*ones(p,1))/dx;     % Frictional loss 
G = (-1)^n * Qvec(1,:) * g * sin(theta);           % Gravitational loss
S = F'+G;                       % Source terms

end
