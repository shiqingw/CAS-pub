function [fplus, fminus] = flux_vec(qvec,A_plus,A_minus)
% Function for flux calculation using vector form

% Get the required parameters from matrix A
% 
% Calculation of Qi & Qi+1 for Fi+1/2
Qip = qvec(:,2:end-1);          % Corresponding Qi
Qip1 = qvec(:,3:end);           % Corresponding Qi+1

% Calculation of Qi & Qi-1 for Fi-1/2
Qim1 = qvec(:,1:end-2);         % Corresponding Qi-1
Qim = qvec(:,2:end-1);          % Corresponding Qi

fplus = A_plus*Qip + A_minus*Qip1;  % Fi+1/2
fminus = A_plus*Qim1 + A_minus*Qim; % Fi-1/2

end
