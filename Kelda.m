function [ F ] = Kelda(q_pump,d,area,dL,rho,mu)
%This function provides the friction function from Kelda

[F, ~] = annularfriction(q_pump,d,area,dL,rho,mu); % Returns pressure drop + dpdq
end

