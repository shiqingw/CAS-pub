function [A_plus, A_minus, L_eigv] = prepare_A(A)

% Prepares different components of matrix A for flux calculation

[R, lambda] = eig(A);                   % Gives right eigen vector R and
                                      	% diagonal eigen values lambda
m = size(lambda);
L_plus = zeros(m);
L_minus = zeros(m);

 for i=1:m(1)
     for j=1:m(2)
         L_plus(i,j) = lambda(i,j)*(lambda(i,j)>0);
         L_minus(i,j) = lambda(i,j)*(lambda(i,j)<0);
     end
 end
                                       
A_plus = R*L_plus/R;       % A plus for flux calculation
A_minus = R*L_minus/R;     % A minus for flux calculation
L_eigv = inv(R);                        % Left EV is inverse of right EV

end

