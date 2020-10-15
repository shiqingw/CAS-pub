function x = newton_solver(F,J,x0,max_itr,tol)
% This function compute the solution of F(x)=0 by using newton method
% arguments of newton_solver:       % F: the function to be solved
                                    % J: Jacobian of F
                                    % x0: the initial condition
                                    % n_it: number of itrations                                
                                    % tol: acceptable tlorance
                                    
x=x0;                               %initial value is assigned 
i=0;                                %counter
res = 1;                            %initial value of residual


while (i<max_itr)&&(res>tol)       
    i=i+1;                            
    x=x-J(x)\F(x);                  %x(i+1)=x(i)+J(x(i))^-1*F(x(1))
    res = norm(F(x),2);             %second norm of F(x(i)
end


end